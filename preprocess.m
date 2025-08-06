% PREPROCESS prepares primer extension sequencing data by extracting
% product and template regions for data meeting quality thresholds.
%
% preprocess(FQ1,FQ2,INDEX,OUT) generates results using the default options 
% for forward read file FQ1, reverse read file FQ2, and index read file 
% INDEX, specified as filenames (strings). The associated files must be 
% sorted and have a 1:1 correspondence of sequencing reads. All outputs are
% generated in the folder OUT, which is created if it does not exist.
%
% SUMMARY=preprocess(...) returns a structure SUMMARY with summary results.
%
% [SUMMARY,OPTIONS]=preprocess(...) also returns a structure OPTIONS that
% includes the full set of options (defaults as well as custom options 
% specified as described below) used in preprocessing.
%
% preprocess(...,'property1',value1,'property2',value2,...) generates 
% results using the additional specified options. Valid options are:
%
% Property          Default Value                    Description
% -------------------------------------------------------------------------
%
% SAMPLE
%
% 'logfile'         fullfile(OUT,'PreProcessed.txt') log file
% 'matfile'         fullfile(OUT,'PreProcessed.mat') results file
% 'map1'            [FQ1 '.map.mat']                 [map files enable
% 'map2'            [FQ2 '.map.mat']                 rapid random access,
% 'index_map        [INDEX '.map.mat']               produced by fastqmap]
%
% PROCESSING
%
% 'reads_per_block' 10000                            files per thread
% 'q_index'         26                               minimum index quality
% 'q'               30                               minimum read quality
% 'q_encoding'      33                               ASCII offset for q = 0
%
% CONSTRUCT
%
% 'fix3'            'GATCGTCGGACTGTAGAACTCTG'        construct fix3 region 
% 'fix2'            'AGATCGGAAGAGCACACGTCTGA'        construct fix2 region
% 'fix1'            'GCATGCGACTAAACGTCGCATGC'        construct fix1 region
% 'template_length' 6                                template region length
% 'prefix'          'TCTA'                           prefix region sequence
%
% Christopher E. Carr (chrisc@mit.edu)
% 
function [summary,options] = preprocess(fq1,fq2,index,out,varargin)
    %% OPTIONS
    % Set any user specified options
    useroptions = args2options(varargin);
    % Set all options to defaults or user specified options
    options = []; % initial empty options
    options = fieldcheck(options,'logfile', fullfile(out,'PreProcessed.txt'), useroptions);
    options = fieldcheck(options,'matfile', fullfile(out,'PreProcessed.mat'), useroptions);
    options = fieldcheck(options,'map1', [fq1 '.map.mat'], useroptions);
    options = fieldcheck(options,'map2', [fq2 '.map.mat'], useroptions);
    options = fieldcheck(options,'reads_per_block', 10000, useroptions);
    % preprocess_read options
    options = fieldcheck(options,'q',             30, useroptions);
    options = fieldcheck(options,'qencoding',     33, useroptions);
    options = fieldcheck(options,'fix3','GATCGTCGGACTGTAGAACTCTG',useroptions);
    options = fieldcheck(options,'fix2','AGATCGGAAGAGCACACGTCTGA',useroptions);
    options = fieldcheck(options,'fix1','GCATGCGACTAAACGTCGCATGC', useroptions);
    options = fieldcheck(options,'template_length', 6, useroptions);
    options = fieldcheck(options,'prefix','TCTA', useroptions);
    options = fieldcheck(options,'index_map',[index '.map.mat'], useroptions);
    options = fieldcheck(options,'q_index',26,useroptions);
    
    processing_options = {'q',options.q,'qencoding',options.qencoding,'fix3',options.fix3,'fix2',options.fix2,'fix1',options.fix1,'template_length',options.template_length,'prefix',options.prefix,'q_index',options.q_index};
    incorporation_bin_edges = [0:1:(options.template_length+2)]-0.5;
    
    %% Verify output folder    
    if ~exist(out,'dir')
        mkdir(out); 
        
        %% Verify mandatory inputs and associated files
        if ~exist(fq1,'file'), error('Input file %s not found',fq1); end
        if ~exist(fq2,'file'), error('Input file %s not found',fq2); end
        if ~exist(index,'file'), error('Input file %s not found',index); end
        if ~exist(options.map1,'file'), error('Map file %s not found',options.map1); end
        if ~exist(options.map2,'file'), error('Map file %s not found',options.map2); end
        if ~exist(options.index_map,'file'), error('Map file %s not found',options.index_map); end
            
        %% Start Logging
        tic; cpu0 = cputime;
        results.options = options;
        fid_log = fopen(options.logfile,'w');
        fprintf(fid_log,'%s\tpreprocess: Started logging\n',datestr(now));
        
        %% Load map files
        map1 = load(options.map1); map1 = map1.map;
        map2 = load(options.map2); map2 = map2.map;
        index_map = load(options.index_map); index_map = index_map.map;
        
        %% B.1.1 Read numbers match.
        % Check number of reads is identical
        if ~(size(map1,1)==size(map2,1)) 
            error('Number of reads in map files %s and %s do not match.',...
                options.map1,options.map2);
        end
        if ~(size(map1,1)==size(index_map,1))
            error('Number of reads in map files %s and %s do not match.',...
            options.map1,options.index_map);
        end        
        % Get number of read pairs
        N_read_pairs = size(map1,1);
        
        %% B.1.2 Read blocks.
        % Calculate number of blocks
        N_blocks = ceil(N_read_pairs/options.reads_per_block);

        parfor block=1:N_blocks
        %for block=1:N_blocks
            %% Get Ready to Process this block
            % Get start for this block
            block_start = (block-1)*options.reads_per_block+1;
            % Get end for this block
            block_end = min(N_read_pairs,block*options.reads_per_block);
            % Create filename for this block
            fn_block_base = fullfile(out,sprintf('block_%03d',block));
            fn_block_mat = [fn_block_base '.mat'];
            fn_block_fasta = [fn_block_base '.fasta'];
            % Extract block files
            [hdr1_b,seq1_b,qual1_b]=ffastqread(fq1,map1(block_start:block_end,:));
            [hdr2_b,seq2_b,qual2_b]=ffastqread(fq2,map2(block_start:block_end,:));
            [hdr3_b,seq3_b,qual3_b]=ffastqread(index,index_map(block_start:block_end,:));
            % Process this set of reads sequentially
            N_reads_in_block = block_end - block_start + 1;
            % preallocate memory for this block
            block_results = struct('bHeadersMatch',NaN,'bIndexQualityPass',NaN,'bR1QualityPass', NaN,...
                'bR2QualityPass', NaN,'bQualityPass', NaN,...
                'bR1Fix2Match', NaN,'bR1Fix1Match', NaN,...
                'bR2Fix1Match', NaN,'bR1PrefixMatch', NaN,'bR1PrefixStart',NaN,...
                'bR1TemplateLengthMatch',NaN,'bR1PrefixLengthMatch',NaN,...
                'bR2TemplateLengthMatch',NaN,'bR2PrefixLengthMatch',NaN,...
                'bR2PrefixMatch', NaN,'bR2Fix3Match',NaN,'bProductLengthPass', NaN,...
                'bProductLengthsMatch', NaN,'bSequencesMatch', NaN,...
                'bDiscard', NaN,'N_incorporations', NaN,...
                'Product','', 'Template', '',...
                'bIncorporation', repmat({NaN},[N_reads_in_block 1]));

            %% Process each read in the block
            for r=1:N_reads_in_block
                H1 = hdr1_b{r}; R1 = seq1_b{r}; Q1 = qual1_b{r};
                H2 = hdr2_b{r}; R2 = seq2_b{r}; Q2 = qual2_b{r};
                H3 = hdr3_b{r}; R3 = seq3_b{r}; Q3 = qual3_b{r};
                % Process this read
                block_results(r) = preprocess_read(H1,R1,Q1,H2,R2,Q2,H3,R3,Q3,processing_options{:});
            end

            %% Calculate block statistics
            block_summary = struct(...
            'ForwardReadFile',fq1,...
            'ReverseReadFile',fq2,...
            'ForwardMapFile',options.map1,...
            'ReverseMapFile',options.map2,...
            'Block',block,...
            'Block_Read_ID_Start',block_start,...
            'Block_Read_ID_Stop',block_end,...
            'Reads',N_reads_in_block,...
            'HeadersMatch',nansum([block_results.bHeadersMatch]),...
            'IndexQualityPass',nansum([block_results.bIndexQualityPass]),...
            'R1QualityPass',nansum([block_results.bR1QualityPass]),...
            'R2QualityPass',nansum([block_results.bR2QualityPass]),...
            'QualityPass',nansum([block_results.bQualityPass]),...
            'R1Fix2Match',nansum([block_results.bR1Fix2Match]),...
            'R1Fix1Match',nansum([block_results.bR1Fix1Match]),...
            'R2Fix1Match',nansum([block_results.bR2Fix1Match]),... 
            'R1TemplateLengthMatch',nansum([block_results.bR1TemplateLengthMatch]),...
            'R1PrefixLengthMatch',nansum([block_results.bR1PrefixLengthMatch]),...
            'R2TemplateLengthMatch',nansum([block_results.bR2TemplateLengthMatch]),...
            'R2PrefixLengthMatch',nansum([block_results.bR2PrefixLengthMatch]),...
            'R1PrefixMatch',nansum([block_results.bR1PrefixMatch]),...
            'R1PrefixStart',nansum([block_results.bR1PrefixStart]),...
            'R2PrefixMatch',nansum([block_results.bR2PrefixMatch]),...
            'R2Fix3Match',nansum([block_results.bR2Fix3Match]),...
            'ProductLengthPass',nansum([block_results.bProductLengthPass]),...
            'ProductLengthsMatch',nansum([block_results.bProductLengthsMatch]),...
            'SequencesMatch',nansum([block_results.bSequencesMatch]),...
            'Discard',nansum([block_results.bDiscard]),...
            'Keep',nansum([block_results.bSequencesMatch]),...
            'ProductNumber',nansum([block_results.bIncorporation]),...
            'ProductLengthHistogram',histcounts([block_results.N_incorporations],incorporation_bin_edges));
            % Save the block results using a helper function to avoid transparency issues inside a parfor loop
            % https://www.mathworks.com/matlabcentral/answers/135285-how-do-i-use-save-with-a-parfor-loop-using-parallel-computing-toolbox
            save_block_results(fn_block_mat,block_summary,block_results);

            %% Save block fasta file
            % Get IDs
            block_read_ids = block_start:block_end;
            % Write to file
            fid_fasta = fopen(fn_block_fasta,'w');
            for r=1:N_reads_in_block
                % Write out results for all valid reads in the block
                if ~isnan(block_results(r).bIncorporation)
                    % Write out fasta as ID with 5' and 3' sequences on
                    % separate lines
                    fprintf(fid_fasta,'>%d\n',block_read_ids(r));
                    fprintf(fid_fasta,'%s\n',block_results(r).Product);
                    fprintf(fid_fasta,'%s\n',block_results(r).Template);
                end
            end
            fclose(fid_fasta);
            %% Log block results to log file
            if isempty(getCurrentTask())
                % only log if we are not running in parallel because if we
                % are running in parallel, file id is invalid
                fprintf(fid_log,'%s\tProcessed block %d (%d reads)\n',datestr(now),block,N_reads_in_block);
            end
        end
        % Summarize statistics across blocks
        for block=1:N_blocks
            fn_block_base = fullfile(out,sprintf('block_%03d',block));
            fn_block_mat = [fn_block_base '.mat'];
            load(fn_block_mat);
            if (block==1)
                summary.ForwardReadFile = fq1;
                summary.ReverseReadFile = fq2;
                summary.IndexReadFile = index;
                summary.ForwardMapFile = options.map1;
                summary.ReverseMapFile = options.map2;
                summary.IndexMapFile = options.index_map;
                summary.Blocks = N_blocks;
                summary.Read_Pairs = N_read_pairs;
                summary.HeadersMatch = 0;
                summary.IndexQualityPass = 0;
                summary.R1QualityPass = 0;
                summary.R2QualityPass = 0;
                summary.QualityPass = 0;
                summary.R1Fix2Match = 0;
                summary.R1Fix1Match = 0;
                summary.R2Fix1Match = 0;
                summary.R1TemplateLengthMatch = 0;
                summary.R1PrefixLengthMatch = 0;
                summary.R2TemplateLengthMatch = 0;
                summary.R2PrefixLengthMatch = 0;
                summary.R1PrefixMatch = 0;
                summary.R1PrefixStart = 0;
                summary.R2PrefixMatch = 0;
                summary.R2Fix3Match = 0;
                summary.ProductLengthPass = 0;
                summary.ProductLengthsMatch = 0;
                summary.SequencesMatch = 0;
                summary.Discard = 0;
                summary.Keep = 0;
                summary.ProductNumber = 0;
                summary.ProductLengthHistogram = zeros(size(incorporation_bin_edges)-[0 1]);
            end

            summary.HeadersMatch = summary.HeadersMatch + block_summary.HeadersMatch;
            summary.IndexQualityPass = summary.IndexQualityPass + block_summary.IndexQualityPass;
            summary.R1QualityPass = summary.R1QualityPass + block_summary.R1QualityPass;
            summary.R2QualityPass = summary.R2QualityPass + block_summary.R2QualityPass;
            summary.QualityPass = summary.QualityPass + block_summary.QualityPass;
            summary.R1Fix2Match = summary.R1Fix2Match + block_summary.R1Fix2Match;
            summary.R1Fix1Match = summary.R1Fix1Match + block_summary.R1Fix1Match;
            summary.R2Fix1Match = summary.R2Fix1Match + block_summary.R2Fix1Match;
            summary.R1TemplateLengthMatch = summary.R1TemplateLengthMatch + block_summary.R1TemplateLengthMatch;
            summary.R1PrefixLengthMatch = summary.R1PrefixLengthMatch + block_summary.R1PrefixLengthMatch;
            summary.R2TemplateLengthMatch = summary.R2TemplateLengthMatch + block_summary.R2TemplateLengthMatch;
            summary.R2PrefixLengthMatch = summary.R2PrefixLengthMatch + block_summary.R2PrefixLengthMatch;
            summary.R1PrefixMatch = summary.R1PrefixMatch + block_summary.R1PrefixMatch;
            summary.R1PrefixStart = summary.R1PrefixStart + block_summary.R1PrefixStart;
            summary.R2PrefixMatch = summary.R2PrefixMatch + block_summary.R2PrefixMatch;
            summary.R2Fix3Match = summary.R2Fix3Match + block_summary.R2Fix3Match;
            summary.ProductLengthPass = summary.ProductLengthPass + block_summary.ProductLengthPass;
            summary.ProductLengthsMatch = summary.ProductLengthsMatch + block_summary.ProductLengthsMatch;
            summary.SequencesMatch = summary.SequencesMatch + block_summary.SequencesMatch;
            summary.Discard = summary.Discard + block_summary.Discard;
            summary.Keep = summary.Keep + block_summary.Keep;
            summary.ProductNumber = summary.ProductNumber + block_summary.ProductNumber;
            summary.ProductLengthHistogram = summary.ProductLengthHistogram + block_summary.ProductLengthHistogram;
        end

        % Save summary of all blocks
        save(options.matfile,'summary','options');        

        %% Make zip archive of all block MAT files
        zip(fullfile(out,'blocks_matfiles.zip'),fullfile(out,'block_*.mat'));
        % remove individual block MAT files
        delete(fullfile(out,'block_*.mat'));

        %% Concatenate all block fasta files
        system(sprintf('cat %s > %s',fullfile(out,'block_*.fasta'),fullfile(out,'PreProcessed.fasta')));
        % remove individual block fasta files
        delete(fullfile(out,'block_*.fasta'));

        %% Calculate execution time
        cpu1 = cputime;
        dt = toc; dcpu = cpu1-cpu0;

        %% Log summary results
        fprintf(fid_log,'%s\tCompleted analysis\n',datestr(now));
        fprintf(fid_log,'\n');
        fprintf(fid_log,'Total Execution Time:\t%g seconds\n',double(dt));
        fprintf(fid_log,'Total CPU Time:      \t%g seconds\n',double(dcpu));
        fprintf(fid_log,'\n');
        fprintf(fid_log,'SUMMARY:\n');
        fprintf(fid_log,'\n');
        fprintf(fid_log,'Dataset:\n');
        fprintf(fid_log,'Forward Read File:\t%s\n',summary.ForwardReadFile);
        fprintf(fid_log,'Reverse Read File:\t%s\n',summary.ReverseReadFile);
        fprintf(fid_log,'Index Read File:  \t%s\n',summary.IndexReadFile);
        fprintf(fid_log,'Forward Map File:\t%s\n',summary.ForwardMapFile);
        fprintf(fid_log,'Reverse Map File:\t%s\n',summary.ReverseMapFile);
        fprintf(fid_log,'Index Map File:  \t%s\n',summary.IndexMapFile);
        fprintf(fid_log,'Read Pairs:\t%d\n',summary.Read_Pairs);
        fprintf(fid_log,'\n');
        fprintf(fid_log,'Options:\n');
        fprintf(fid_log,'Minimum Quality (Q):\t%d\n',options.q);
        fprintf(fid_log,'Minimum Index Quality (Q):\t%d\n',options.q_index);
        fprintf(fid_log,'Quality Encoding:\t%d\n',options.qencoding);
        fprintf(fid_log,'Prefix Sequence:\t%s\n',options.prefix);
        fprintf(fid_log,'Template Length:\t%d\n',options.template_length);
        fprintf(fid_log,'Fix1 Sequence:\t%s\n',options.fix1);
        fprintf(fid_log,'Fix2 Sequence:\t%s\n',options.fix2);
        fprintf(fid_log,'Fix3 Sequence:\t%s\n',options.fix3);
        fprintf(fid_log,'\n');
        fprintf(fid_log,'Results:\n');
        fprintf(fid_log,'Headers Match:\t%d\n',summary.HeadersMatch);
        fprintf(fid_log,'Index Quality Match: \t%d\n',summary.IndexQualityPass);
        fprintf(fid_log,'Read 1 Quality Match:\t%d\n',summary.R1QualityPass);
        fprintf(fid_log,'Read 2 Quality Match:\t%d\n',summary.R2QualityPass);
        fprintf(fid_log,'Quality Pass:\t%d\n',summary.QualityPass);
        fprintf(fid_log,'Read 1 Fix 2 Match:\t%d\n',summary.R1Fix2Match);
        fprintf(fid_log,'Read 1 Fix 1 Match:\t%d\n',summary.R1Fix1Match);
        fprintf(fid_log,'Read 2 Fix 1 Match:\t%d\n',summary.R2Fix1Match);
        fprintf(fid_log,'Read 1 Template Length Match:\t%d\n',summary.R1TemplateLengthMatch);
        fprintf(fid_log,'Read 1 Prefix Length Match:\t%d\n',summary.R1PrefixLengthMatch);
        fprintf(fid_log,'Read 1 Prefix Start:\t%d\n',summary.R1PrefixStart);
        fprintf(fid_log,'Read 2 Template Length Match:\t%d\n',summary.R2TemplateLengthMatch);
        fprintf(fid_log,'Read 2 Prefix Length Match:\t%d\n',summary.R2PrefixLengthMatch);
        fprintf(fid_log,'Read 1 Prefix Match:\t%d\n',summary.R1PrefixMatch);
        fprintf(fid_log,'Read 2 Prefix Match:\t%d\n',summary.R2PrefixMatch);
        fprintf(fid_log,'Read 2 Fix3 Match:\t%d\n',summary.R2Fix3Match);
        fprintf(fid_log,'Product Length Pass:\t%d\n',summary.ProductLengthPass);
        fprintf(fid_log,'Product Lengths Match:\t%d\n',summary.ProductLengthsMatch);
        fprintf(fid_log,'Read 1/Read 2 Sequences Match:\t%d\n',summary.SequencesMatch);
        fprintf(fid_log,'Read Pair Discarded:\t%d\n',summary.Discard);
        fprintf(fid_log,'Read Pair Retained:\t%d\n',summary.Keep);
        fprintf(fid_log,'Read Pairs with >=1 Incorporations:\t%d\n',summary.ProductNumber);
        fprintf(fid_log,'Product Length Histogram (0-n):\t[');
        inc_hist_tmp = sprintf('%d ',summary.ProductLengthHistogram);
        fprintf(fid_log,'%s]\n',inc_hist_tmp(1:end-1));
        fclose(fid_log);
    else
        % Output folder exists
        fprintf('Output folder %s exists.\n',out);
        fprintf('Please remove existing folder or modify output folder name and rerun.\n\n',out);
    end
end

