% CHARACTERIZE analyzes matched product-template regions of primer 
% extension sequencing data extracted via preprocessing with PREPROCESS.
%
% characterize(OUT) stores the results in the folder OUT.
%
% S=characterize(OUT) returns a structure S with summary 
% results.
%
% [S,OP]=characterize(OUT) also returns a structure OP that includes the
% full set of options (defaults as well as custom options specified as
% described below) used in characterization.
%
% [S,OP,POP]=characterize(OUT) also returns a structure POP that the
% options used for preprocessing.
%
% characterize(OUT,'property1',value1,'property2',value2,...) generates 
% results using the additional specified options. Valid options are:
%
% Property          Default Value                    Description
% -------------------------------------------------------------------------
%
% SAMPLE
%
% 'logfile'         fullfile(OUT,'Characterized.txt')  log file
% 'matfile'         fullfile(OUT,'PreProcessed.mat')   preprocessing MAT
% 'fasta'           fullfile(OUT,'PreProcessed.fasta') preprocessing fasta
% 'save_matfile'    fullfile(OUT,'Characterized.mat')  characterize results
% 'plots'           fullfile(OUT,'plots')              folder for plots
% 'extra'           fullfile(OUT,'extra')              folder for extra files
% 'control'                                            see below
% 'normalization_factor'                               see below
%
% 'control': MAT filename of any associated control condition. For internal
% normalization, use ''. For normalization with a control condition, use
% output of running characterize() on this control condition, e.g. MAT file
% such as './out/Control/Characterized.mat'.
%
% 'normalization_factor' The normalization factor is a matrix of base 
% frequencies with four rows, corresponding to A, C, G, U, and one column 
% for each template position. This can be specified as a matrix of values
% explicitely, or as a matlab expression, such as:
% [[0.25*ones(4,6)],[1;0;0;0]]
%
% This can also be specified as a sequence, and ambiguous nucleotides can
% be specified. For example, the above expression is equivalent to:
% 'NNNNNNA'
%
% This example reflects a template region synthesized to be a uniform
% random distribution across bases {A, C, G, U} and with a defined terminal
% A base just before caged base(s). Formally, the template length is the
% six N bases, and the A is the first base of the construct outside the
% template region.
%
% The normalization factor is used to adjust base frequencies to compensate
% for synthesizer errors or selection bias in NERPE experiments.
%
% A template normalization factor is computed as the ratio of the 
% normalization_factor to the template frequencies from either the sample
% itself or from a control, if a control file is specified using 'control'.
%
% In the absence of synthesizer errors or selection bias, the resulting
% template normalization factor will be unity. Deviations from unity thus
% measure synthesizer error and/or selection bias. 
%
% Base frequencies are normalized by this template normalization factor to 
% correct for these biases.
%
% Christopher E. Carr (chrisc@mit.edu)
% 

function [summary,options,preprocessing_options] = characterize(out,varargin)
    %% OPTIONS
    % Set any user specified options
    useroptions = args2options(varargin);
    % Set all options to defaults or user specified options
    options = []; % initial empty options
    options = fieldcheck(options,'logfile', fullfile(out,'Characterized.txt'), useroptions);
    options = fieldcheck(options,'matfile', fullfile(out,'PreProcessed.mat'), useroptions);
    options = fieldcheck(options,'fasta', fullfile(out,'PreProcessed.fasta'), useroptions);
    options = fieldcheck(options,'save_matfile', fullfile(out,'Characterized.mat'), useroptions);
    options = fieldcheck(options,'plots',fullfile(out,'plots'), useroptions);
    options = fieldcheck(options,'extra',fullfile(out,'extra'), useroptions);
    options = fieldcheck(options,'control','', useroptions);
    options = fieldcheck(options,'internal_normalization',0, useroptions);
    
    % preprocess_read options are obtained from MAT file
    options = fieldcheck(options,'normalization_factor',[[0.25*ones(4,6)],[1;0;0;0]], useroptions);
    
    % Handle NaN value provided as options.control
    if ~ischar(options.control), options.control = ''; end
    
    % check for valid normalization
    if strcmp(options.control,'')
        if options.internal_normalization==0
            % Warn user, using internal normalization because control was not specified
            fprintf('characterize: using internal normalization because no control was specified (this is a normal outcome when processing a control)\n')
            options.internal_normalization = 1;
        end
    end

    if ~exist(out,'dir')
        error('characterize: Preprocessing output directory %s not found.',out);
    else
        % Found output directory
        
        %% Verify mandatory inputs and associated files
        if ~exist(options.matfile,'file'), error('Preprocessing MAT file %s not found',options.matfile); end
        if ~exist(options.fasta,'file'), error('Preprocessed Fasta file %s not found',options.fasta); end
        if ~isempty(options.control), if ~exist(options.control,'file')
                error('Conrol MAT file %s not found',options.control);
            end
        end
        %% Start Logging
        tic; cpu0 = cputime;
        fid_log = fopen(options.logfile,'a');
        fprintf(fid_log,'%s\tcharacterize: Started logging\n',datestr(now));
    
        % Load the MATFILE (already verified to exist)
        matfile = load(options.matfile);
        if ~isempty(options.control)
            control = load(options.control);
        end
        % Initialize other outputs
        summary = matfile.summary;
        preprocessing_options = matfile.options;
        
        % Define way to map RNA ASCII to RNA base
        ascii_2_baseid = sparse([1 1 1 1 1],[65 67 71 85 45],[1 2 3 4 5]);
        % Define columns for position in template or product
        col_id = 1:(matfile.options.template_length+1);

        % Load the preprocessed fasta file; its not a standard fasta format
        % so we need to call the helper routine.
        [H,P,T]=read_nerpe_fasta(options.fasta);
        
        % Convert * to - in product (- is standard for deletion)
        P = strrep(P,'*','-');
                
        %% C.1 Extension Events
        
        % C.1.1 Calculate the raw number of unextended, +1 extended, 
        % +2 extended, etc. blocks. This includes everything--
        % complementary, mismatch, unextended--and is meant to be akin to 
        % what we would learn by looking at a gel.
        
        % This was previously computed as summary.IncorporationHistogram
        % but we're doing it again in case we are dealing with filtered
        % fasta file.
        
        % Get appropriate bin edges
        incorporation_bin_edges = [0:1:(preprocessing_options.template_length+2)]-0.5;
        
        % Initialize
        summary.Incorporation = 0;
        summary.IncorporationHistogram = zeros(size(incorporation_bin_edges)-[0 1]);
        summary.IncorporationHistogram_Relative = zeros(size(incorporation_bin_edges)-[0 1]);
        
        % Calculate
        p_len = cellfun(@numel,strrep(P,'-',''));
        summary.Incorporation = sum(p_len>0);
        summary.IncorporationHistogram = histcounts(p_len,incorporation_bin_edges);
        
        % C.1.2 Relative-to-total unextended/extension events
        summary.IncorporationHistogram_Relative = ...
            summary.IncorporationHistogram ./ ...
            sum(summary.IncorporationHistogram);
        
        %% C.2 Template Composition 
        
        % C.2.1 Raw raw position-dependent base composition of template
        T_Counts = GetBaseCounts(T,5,ascii_2_baseid,matfile.options.template_length+1);
        % Save to summary
        summary.TemplateCounts = T_Counts;
        
        % C.2.2 Position-dependent base composition of the template as
        % frequencies
        T_Freq = T_Counts(1:4,:)./repmat(nansum(T_Counts(1:4,:)),[4 1]);
        % Save to summary
        summary.TemplateFrequencies = T_Freq;
        
        % C.2.3 Compare the position-dependent base composition of the 
        % template to that of any specified control template (if 
        % specified). Simply compare one to the other by taking the 
        % ratios of the bases at each position (raw or frequencies). 
        % (e.g.: T1 in experimental = 20% G, T1 in control = 23% G; 
        % .2/.23 = 0.87.) The null hypothesis is that the ratios will 
        % equal 1, whereas any deviation from unity could point to a 
        % sequence-dependent 2'-5' linkage signature.
        if ~isempty(options.control)
            T_Freq_vs_Control = T_Freq./control.summary.TemplateFrequencies;
            summary.TemplateFrequenciesVersusControl = T_Freq_vs_Control;
        else
            summary.TemplateFrequenciesVersusControl = NaN(size(T_Freq));
        end
        
        % C.2.4 Estimate normalization factors to compensate for variation
        % in template composition.
        
        % Always use an external control template for normalization unless 
        % internal_normalization is set to nonzero
        
        if options.internal_normalization
            % Internal
            T_normalization_factors = options.normalization_factor ./ T_Freq;
        else
            % External (control)
            T_normalization_factors = options.normalization_factor ./ control.summary.TemplateFrequencies;
        end
        T_normalization_factors(isinf(T_normalization_factors))=NaN;

        % Save to summary
        summary.TemplateNormalizationFactors = T_normalization_factors;
        
        %% C.3 Product Grouping
        
        % P_Composition = GetBaseCounts(P,5,ascii_2_baseid,matfile.options.template_length+1);
        
        % C.3.1 Group Complementary, Mismatch, and Unextended sets
        
        % preallocate memory
        bCompSet = logical(size(P));
        bMMSet = logical(size(P));
        bUnSet = logical(size(P));
        % parallel execution
        parfor k=1:max(size(P))
        % serial execution
        %for k=1:max(size(P))
            P_k = strrep(P{k},'-','');
            n_k = numel(P_k);
            if n_k==0
                bCompSet(k) = false;
                bMMSet(k) = false;
                bUnSet(k) = true;
            else
                % Either Mismatch or Complementary Set
                % Note: template is stored as 3'-5' so we don't need to take
                % the reverse complement, just the complement.
                T_k = T{k}; T_k = dna2rna(seqcomplement(rna2dna(T_k(1:n_k))));
                if strcmp(P_k,T_k)
                    bCompSet(k) = true;
                    bMMSet(k) = false;
                    bUnSet(k) = false; 
                else
                    bCompSet(k) = false;
                    bMMSet(k) = true;
                    bUnSet(k) = false;
                end
            end
        end
        
        % Get Comp Set and MM Set for Products and Templates
        P_CompSet = P(bCompSet);
        P_MMSet = P(bMMSet);
        P_UnSet = P(bUnSet);
        T_CompSet = T(bCompSet);
        T_MMSet = T(bMMSet);
        T_UnSet = T(bUnSet);
       
        % C.3.2 Count the number of blocks (alignments) in each set
        N_P_CompSet = sum(bCompSet);
        N_P_MMSet = sum(bMMSet);
        N_P_UnSet = sum(bUnSet);
        
        bBlockCountsMatch = (N_P_CompSet+N_P_MMSet+N_P_UnSet)==numel(P);
        
        % Save results to summary
        summary.bBlockCountsMatch = bBlockCountsMatch;
        summary.ProductComplementarySetN = N_P_CompSet;
        summary.ProductMismatchSetN = N_P_MMSet;
        summary.ProductUnincorporatedSetN = N_P_UnSet;
        
        %% Product features in the Complementary Set and Unextended Set
        
        % C.4.1 Position dependent raw base composition of product
        
        % Complementary set
        P_CompSet_Counts = GetBaseCounts(P_CompSet,5,ascii_2_baseid,...
            matfile.options.template_length+1);
        % Save to summary
        summary.ProductComplementarySetCounts = P_CompSet_Counts;
        
        % Unextended Set
        P_UnSet_Counts = GetBaseCounts(P_UnSet,5,ascii_2_baseid,...
            matfile.options.template_length+1);
        % Save to summary
        summary.ProductUnincorporatedSetCounts = P_UnSet_Counts;
        
        % C.4.2 Complementary Set. Normalize the raw position-dependent 
        % base composition of the product (C.4.1) to the template using 
        % the normalization factors (C.2.4), excluding nulls.
        
        % To normalize, need to permute rows of template normalization
        % factor computed in step C.2.4
        %
        %  Product   Template           Row of template
        %  A         U                  4
        %  C         G                  3
        %  G         C                  2
        %  U         A                  1
        %  -         (not possible)     (not applicable)
        %
        
        % Normalization factor
        P_norm_factor = T_normalization_factors([4 3 2 1],:);
        
        % Compute normalized composition for A, C, G, and U
        P_CompSet_Counts_Normalized = P_CompSet_Counts(1:4,:) .* P_norm_factor;
        % Save to summary
        summary.ProductComplementarySetCountsNormalized = P_CompSet_Counts_Normalized;
        
        % C.4.3 Complementary Set. Compute the normalized position-dependent 
        % base composition of the product (C.4.2) as frequencies.
        P_CompSet_Freq_Normalized = P_CompSet_Counts_Normalized./repmat(nansum(P_CompSet_Counts_Normalized),[4 1]);
        % Save to summary
        summary.ProductComplementarySetFrequenciesNormalized = P_CompSet_Freq_Normalized;
        
        % C.4.4 and C.4.5
        if ~isempty(options.control)
            % C.4.4 Complementary Set. Normalize the raw position-dependent 
            % base composition of the product (C.4.1) to the control template 
            % using the control template normalization factors (C.2.4).
            P_control_factor = control.summary.TemplateNormalizationFactors([4 3 2 1],:);
            P_CompSet_Counts_ControlTemplateNormalized = P_CompSet_Counts(1:4,:) .* P_control_factor;    
            % Save to summary
            summary.ProductComplementarySetCountsControlTemplateNormalized = P_CompSet_Counts_ControlTemplateNormalized;
            % C.4.5 Complementary Set. List the control-template-normalized 
            % position-dependent base composition of the product (C.4.4) 
            % as frequencies.
            P_CompSet_Freq_ControlTemplateNormalized = P_CompSet_Counts_ControlTemplateNormalized./repmat(nansum(P_CompSet_Counts_ControlTemplateNormalized),[4 1]);
            % Save to summary
            summary.ProductComplementarySetFrequenciesControlTemplateNormalized = P_CompSet_Freq_ControlTemplateNormalized;        
        else
            summary.ProductComplementarySetCountsControlTemplateNormalized = NaN(4,size(P_CompSet_Counts,2));
            summary.ProductComplementarySetFrequenciesControlTemplateNormalized = NaN(4,size(P_CompSet_Counts,2));
        end
        
        %% C.5 Sequence-dependence of nulls (PE termination) in the 
        % Complementary Set and Unextended Set. This analysis allows us to 
        % answer the question: Do any terminal product bases or templating 
        % bases tend to stop PE?
        
        % C.5.1 Complementary Set. Compute the raw position-dependent base 
        % distribution of all terminal product bases (i.e.: the base 
        % distribution of the final base of all products at each position)
        n_max = matfile.options.template_length+1;
        P_Terminal_Product_Base = repmat({repmat('-',[1 n_max])},[numel(P_CompSet) 1]);
        for k=1:numel(P_CompSet)
            P_k = strrep(P_CompSet{k},'-','');
            n_k = numel(P_k);
            if and(n_k<n_max,n_k>0)
                % At least one incorporation, and
                % Length is less than maximum length, implying a deletion
                % Determine final product base sequence
                P_k_terminal_product_base = repmat('-',[1 n_max]);
                P_k_terminal_product_base(n_k)=P_k(n_k);
                P_Terminal_Product_Base{k} = P_k_terminal_product_base;
            end
        end
        % Calculate counts for terminal product bases
        P_Terminal_Product_Base_Counts = GetBaseCounts(P_Terminal_Product_Base,5,ascii_2_baseid,...
            matfile.options.template_length+1);
        % Deletion row is not valid
        P_Terminal_Product_Base_Counts = P_Terminal_Product_Base_Counts(1:4,:);
        % Save to summary
        summary.TerminalProductBaseCounts = P_Terminal_Product_Base_Counts;
        
        % C.5.2 Complementary Set. Normalize the raw position-dependent 
        % base distribution of all terminal product bases (C.5.1) to the 
        % template using the normalization factors (C.2.4).
        
        % Product normalization factor was calculated above as P_norm_factor
        P_Terminal_Product_Base_Counts_Normalized = P_Terminal_Product_Base_Counts.*P_norm_factor;
        % Save to summary
        summary.TerminalProductBaseCountsNormalized = P_Terminal_Product_Base_Counts_Normalized;
        
        % C.5.3 Complementary Set. Compute the normalized position-dependent 
        % base distribution of all terminal product bases (C.5.2) as 
        % frequencies.
        P_Terminal_Product_Base_Freq_Normalized = P_Terminal_Product_Base_Counts_Normalized./repmat(nansum(P_Terminal_Product_Base_Counts_Normalized),[4 1]);  
        % Save to summary
        summary.TerminalProductBaseFrequenciesNormalized = P_Terminal_Product_Base_Freq_Normalized;
        
        % C.5.4 Complementary Set and Unextended Set. Compute the raw 
        % position-dependent base distribution of all templating bases 
        % one base downstream of all terminal products (which, formally, 
        % in this case, includes all unextended products).
        P_CompUnSet = [P_CompSet P_UnSet]; T_CompUnSet = [T_CompSet T_UnSet];
        T_First_Null_Base = repmat({repmat('-',[1 n_max])},[numel(P_CompUnSet) 1]);
        n_max = matfile.options.template_length+1;
        for k=1:numel(P_CompUnSet)
            P_k = strrep(P_CompUnSet{k},'-','');
            n_k = numel(P_k); 
            T_k = T_CompUnSet{k};
            if n_k<n_max
                % Length is less than maximum length, implying a deletion
                % Determine final product base sequence
                T_k_first_null_base = repmat('-',[1 n_max]);
                T_k_first_null_base(n_k+1)=T_k(n_k+1);
                T_First_Null_Base{k} = T_k_first_null_base;
            end
        end
        % Calculate counts 
        T_First_Null_Base_Counts = GetBaseCounts(T_First_Null_Base,5,ascii_2_baseid,...
            matfile.options.template_length+1);
        % Deletion row is not valid
        T_First_Null_Base_Counts = T_First_Null_Base_Counts(1:4,:);
        % Save to summary
        summary.TemplateFirstNullProductBaseCounts = T_First_Null_Base_Counts;
        
        % C.5.5 Complementary Set and Unextended Set. Normalize the 
        % raw position-dependent base distribution of all templating bases 
        % one base downstream of all terminal products (C.5.4) to the 
        % template using the normalization factors (C.2.4).
        T_First_Null_Base_Counts_Normalized = T_First_Null_Base_Counts .* T_normalization_factors;
        % Save to summary
        summary.TemplateFirstNullProductBaseCountsNormalized = T_First_Null_Base_Counts_Normalized;
        
        % C.5.6 Complementary Set and Unextended Set. Compute the 
        % normalized position-dependent base distribution of all 
        % templating bases one base downstream of all terminal 
        % products (C.5.5) as frequencies.
        T_First_Null_Base_Freq_Normalized = T_First_Null_Base_Counts_Normalized./repmat(sum(T_First_Null_Base_Counts_Normalized,1),[4 1]);      
        % Save to summary
        summary.TemplateFirstNullProductFrequenciesNormalized = T_First_Null_Base_Freq_Normalized;
        
        % C.5.7 Complementary Set. Compute the number of template dimers as 
        % a function of incorporation position i, where each dimer is the 
        % template sequence at positions [i i+1].
        
        dimers = {'AA' 'AC' 'AG' 'AU' 'CA' 'CC' 'CG' 'CU' 'GA' 'GC' 'GG' 'GU' 'UA' 'UC' 'UG' 'UU'};
        n_cols = preprocessing_options.template_length;
        TDSeqSpace = NaN(numel(dimers),n_cols);
        for k=1:n_cols
            % Look at position k
            T_k = cellfun(@(x)(x(k)),T_CompSet)';
            T_kp1 = cellfun(@(x)(x(k+1)),T_CompSet)';
            P_k = cellfun(@(x)(x(k)),P_CompSet)';
            bProduct = ~(P_k=='-');
            % Restrict template dimers to positions where there is a product
            T_k = T_k(bProduct);
            T_kp1 = T_kp1(bProduct);
            % Sum up number of each type of dimer
            Pair = cellstr([T_k T_kp1]);
            for j=1:numel(dimers)
                rc_j = sum(cellfun(@(x)(strcmp(x,dimers{j})),Pair));
                TDSeqSpace(j,k)=rc_j;
            end
        end
        % 
        % Save to summary
        summary.TemplateCompSetDimerCounts = TDSeqSpace;
        
        % Frequencies
        TDSeqSpace_Freq = TDSeqSpace./repmat(sum(TDSeqSpace,1),[16 1]);      
        % Save to summary
        summary.TemplateCompSetDimerFrequencies = TDSeqSpace_Freq;
        
        % C.5.8. Compute the number of template dimers (entire set) as 
        % a function of position i, where each dimer is the 
        % template sequence at positions [i i+1].
        
        dimers = {'AA' 'AC' 'AG' 'AU' 'CA' 'CC' 'CG' 'CU' 'GA' 'GC' 'GG' 'GU' 'UA' 'UC' 'UG' 'UU'};
        n_cols = preprocessing_options.template_length;
        TSeqSpace = NaN(numel(dimers),n_cols);
        for k=1:n_cols
            % Look at position k
            T_k = cellfun(@(x)(x(k)),T)';
            T_kp1 = cellfun(@(x)(x(k+1)),T)';
            % Sum up number of each type of dimer
            Pair = cellstr([T_k T_kp1]);
            for j=1:numel(dimers)
                rc_j = sum(cellfun(@(x)(strcmp(x,dimers{j})),Pair));
                TSeqSpace(j,k)=rc_j;
            end
        end
        % 
        % Save to summary
        summary.TemplateDimerCounts = TSeqSpace;
        
        % Frequencies
        TSeqSpace_Freq = TSeqSpace./repmat(sum(TSeqSpace,1),[16 1]);      
        % Save to summary
        summary.TemplateDimerFrequencies = TSeqSpace_Freq;
        
        % C.6 MM space (MM set). What are the properties of products 
        % among MisMatches?
        
        % C.6.1 MisMatch Set. How many products, relative to the total 
        % of all extended products (Complementary Set + MisMatch Set), 
        % have 1, 2, 3, etc. mismatches?
        
        % C.6.2 MisMatch Set. What is the position-dependent distribution 
        % of a MM being terminal?
        
        bMM = repmat(repmat(NaN,[1 n_max]),[N_P_MMSet 1]);
        bTerminal = repmat(repmat(NaN,[1 n_max]),[N_P_MMSet 1]);
                
        for k=1:N_P_MMSet
            P_MM_k = strrep(P_MMSet{k},'-','');
            n_k = numel(P_MM_k);
            T_MM_k = T_MMSet{k};
            % Note: template is stored as 3'-5' so we don't need to take
            % the reverse complement, just the complement.
            T_MM_k = dna2rna(seqcomplement(rna2dna(T_MM_k(1:n_k))));
            % character by character comparison, 1 or 0
            bMM_k = ~eq(P_MM_k,T_MM_k);
            % Store character by character comparison
            bMM(k,:) = [bMM_k repmat(NaN,[1 n_max-n_k])];
            % is it terminal (is a 1 in the last position)?
            bTerminal(k,n_k) = bMM_k(n_k);
        end
        % sum up within a given read to get number of mismatches
        nMM = nansum(bMM,2);
        pMM = nansum(bMM,1);
        pTerminal = nansum(bTerminal,1);
        % Terminal mismatches are not defined in final position due to no
        % knowledge of whether a deletion follows
        pTerminal(end)=NaN;
        % Save results to summary
        bin_edges = [0:1:(preprocessing_options.template_length+2)]-0.5;
        summary.ProductMismatchSetMismatchCounts = pMM;
        summary.ProductMismatchSetTerminalMismatchCounts = pTerminal;
        summary.ProductMismatchSetTerminalMismatchCountsRelativeToAllMismatches = pTerminal./pMM;
        summary.ProductMismatchSetMismatchCountsPerProductHistogram = histcounts(nMM,bin_edges);
        % Normalize to number of bases that have k or more incorporations
        tmp = flip(cumsum(flip(summary.IncorporationHistogram(2:end))));
        summary.ProductMismatchSetMismatchCountsRelativeToAllIncorporations = pMM./tmp;
        
        % C.6.3 MisMatch Set. What is the position-dependent distribution 
        % of a MM followed by a MM?
        
        % Look for [1 1] subsequence indicating [mismatch mismatch]
        pat = [1 1]; n_pat = numel(pat); 
        n_col = (preprocessing_options.template_length+2-numel(pat));
        mm_mm = NaN(N_P_MMSet,n_col);
        for k=1:n_col
            bMM_k = bMM(:,k:(k+numel(pat)-1));
            for j=1:N_P_MMSet
                mm_mm(j,k) = sum(eq(bMM_k(j,:),pat))==n_pat;
            end
        end
        mm_mm_count = nansum(mm_mm,1);    
        % Save to summary
        summary.ProductMismatchSetMismatchFollowedByMismatchCounts = mm_mm_count;
            
        % C.6.4 MisMatch Set. What is the position-dependent distribution 
        % of a MM followed by a correctly paired base?

        % Look for [1 0] subsequence indicating [mismatch non-mismatch]
        pat = [1 0]; n_pat = numel(pat); 
        n_col = (preprocessing_options.template_length+2-numel(pat));
        mm_c = NaN(N_P_MMSet,n_col);
        for k=1:n_col
            bMM_k = bMM(:,k:(k+numel(pat)-1));
            for j=1:N_P_MMSet
                mm_c(j,k) = sum(eq(bMM_k(j,:),pat))==n_pat;
            end
        end
        mm_c_count = nansum(mm_c,1);    
        % Save to summary
        summary.ProductMismatchSetMismatchFollowedByNonMismatchCounts = mm_c_count;
        
        % C.6.5 MisMatch Set. Sequence space of mismatches. Position-
        % dependent prevalence of each possible mismatch (including 
        % orientation!) at each position.
        
        % Define each possibility as a bin at each position.
        % e.g., [Template-Product]... AA, AC, AG, etc.
        % (note complimentary combinations, e.g., AU, will be zero)
        
        % Rows vs. Columns Position 1 2 3 4 ...
        % AA
        % AC
        % AG
        % CA 
        % CC
        % CU
        % GA
        % GG
        % GU
        % UC
        % UG
        % UU
        
        rows = {'AA' 'AC' 'AG' 'CA' 'CC' 'CU' 'GA' 'GG' 'GU' 'UC' 'UG' 'UU'};
        n_cols = preprocessing_options.template_length+1;
        MMSeqSpace = NaN(numel(rows),n_cols);
        for k=1:n_cols
            % Look at mismatches in position k
            T_k = cellfun(@(x)(x(k)),T_MMSet)';
            P_k = cellfun(@(x)(x(k)),P_MMSet)';
            if bitand(~isempty(T_k),~isempty(P_k))
                Pair = cellstr([T_k P_k]);
                for j=1:numel(rows)
                    rc_j = sum(cellfun(@(x)(strcmp(x,rows{j})),Pair));
                    MMSeqSpace(j,k)=rc_j;
                end
            end
        end
        % Add to summary
        summary.ProductMismatchSetMismatchSequenceSpaceCounts = MMSeqSpace;
        
        % To normalize, apply normalization 
        % factors across entire contents of each base-defined set of 
        % bins (i.e.: In position 1, the "A bin" will contain AA, AG, 
        % and AC in some proportion of raw counts of each. Normalize 
        % all 3 bins to the appropriate position 1 normalization factor 
        % (which in this case is the template normalization factor).
        MMSeqSpace_Norm_Factor = T_normalization_factors([1 1 1 2 2 2 3 3 3 4 4 4],:);
        % Normalize
        MMSeqSpace_Normalized = MMSeqSpace.*MMSeqSpace_Norm_Factor;
        % Save to summary
        summary.ProductMismatchSetMismatchSequenceSpaceCountsNormalized = MMSeqSpace_Normalized;
        
        % Normalized frequencies
        MMSeqSpace_Freq_Normalized = MMSeqSpace_Normalized./repmat(sum(MMSeqSpace_Normalized,1),[12 1]);      
        % Save to summary
        summary.ProductMismatchSetMismatchSequenceSpaceFrequenciesNormalized = MMSeqSpace_Freq_Normalized;
        
        %% Context of a mismatch
        
        % Determine base 0 (p-1 base when position is p)
        base0_product = preprocessing_options.fix1(end);
        
        % get each mismatch with context (mwc) and save to file
        if ~exist(options.extra,'dir'), mkdir(options.extra); end
        mwc_saveas = fullfile(options.extra,'MismatchContext.xlsx');
        mwc = mismatch_context(P_MMSet,T_MMSet,base0_product,mwc_saveas);
        
        %% Save results back to matfile
        save(options.save_matfile,'summary','options','preprocessing_options');        
        
        %% Calculate execution time
        cpu1 = cputime;
        dt = toc; dcpu = cpu1-cpu0;

        %% Log summary results and Stop Logging
        fprintf(fid_log,'%s\tcharacterize: Completed analysis\n',datestr(now));
        fprintf(fid_log,'\n');
        fprintf(fid_log,'Total Execution Time:\t%g seconds\n',double(dt));
        fprintf(fid_log,'Total CPU Time:      \t%g seconds\n',double(dcpu));
        fprintf(fid_log,'\n');
        fprintf(fid_log,'SUMMARY:\n');
        fprintf(fid_log,'\n');
        fprintf(fid_log,'Dataset:\n');
        fprintf(fid_log,'Forward Read File:\t%s\n',matfile.summary.ForwardReadFile);
        fprintf(fid_log,'Reverse Read File:\t%s\n',matfile.summary.ReverseReadFile);
        fprintf(fid_log,'Index Read File:  \t%s\n',matfile.summary.IndexReadFile);
        fprintf(fid_log,'Forward Map File:\t%s\n',matfile.summary.ForwardMapFile);
        fprintf(fid_log,'Reverse Map File:\t%s\n',matfile.summary.ReverseMapFile);
        fprintf(fid_log,'Index Map File:  \t%s\n',matfile.summary.IndexMapFile);
        fprintf(fid_log,'Read Pairs:\t%d\n',matfile.summary.Read_Pairs);
        fprintf(fid_log,'\n');
        fprintf(fid_log,'Options:\n');
        fprintf(fid_log,'Minimum Quality (Q):\t%d\n',matfile.options.q);
        fprintf(fid_log,'Minimum Index Quality (Q):\t%d\n',matfile.options.q_index);
        fprintf(fid_log,'Quality Encoding:\t%d\n',matfile.options.qencoding);
        fprintf(fid_log,'Prefix Sequence:\t%s\n',matfile.options.prefix);
        fprintf(fid_log,'Template Length:\t%d\n',matfile.options.template_length);
        fprintf(fid_log,'Fix1 Sequence:\t%s\n',matfile.options.fix1);
        fprintf(fid_log,'Fix2 Sequence:\t%s\n',matfile.options.fix2);
        fprintf(fid_log,'\n');
        fprintf(fid_log,'Preprocessing Results:\n');
        fprintf(fid_log,'Headers Match:\t%d\n',matfile.summary.HeadersMatch);
        fprintf(fid_log,'Index Quality Match: \t%d\n',matfile.summary.IndexQualityPass);
        fprintf(fid_log,'Read 1 Quality Match:\t%d\n',matfile.summary.R1QualityPass);
        fprintf(fid_log,'Read 2 Quality Match:\t%d\n',matfile.summary.R2QualityPass);
        fprintf(fid_log,'Quality Pass:\t%d\n',matfile.summary.QualityPass);
        fprintf(fid_log,'Read 1 Fix 2 Match:\t%d\n',matfile.summary.R1Fix2Match);
        fprintf(fid_log,'Read 1 Fix 1 Match:\t%d\n',matfile.summary.R1Fix1Match);
        fprintf(fid_log,'Read 2 Fix 1 Match:\t%d\n',matfile.summary.R2Fix1Match);
        fprintf(fid_log,'Read 1 Template Length Match:\t%d\n',matfile.summary.R1TemplateLengthMatch);
        fprintf(fid_log,'Read 1 Prefix Length Match:\t%d\n',matfile.summary.R1PrefixLengthMatch);
        fprintf(fid_log,'Read 2 Template Length Match:\t%d\n',matfile.summary.R2TemplateLengthMatch);
        fprintf(fid_log,'Read 2 Prefix Length Match:\t%d\n',matfile.summary.R2PrefixLengthMatch);
        fprintf(fid_log,'Read 1 Prefix Match:\t%d\n',matfile.summary.R1PrefixMatch);
        fprintf(fid_log,'Read 2 Prefix Match:\t%d\n',matfile.summary.R2PrefixMatch);
        fprintf(fid_log,'Product Length Pass:\t%d\n',matfile.summary.ProductLengthPass);
        fprintf(fid_log,'Product Lengths Match:\t%d\n',matfile.summary.ProductLengthsMatch);
        fprintf(fid_log,'Read 1/Read 2 Sequences Match:\t%d\n',matfile.summary.SequencesMatch);
        fprintf(fid_log,'Read Pair Discarded:\t%d\n',matfile.summary.Discard);
        fprintf(fid_log,'Read Pair Retained:\t%d\n',matfile.summary.Keep);        
        fprintf(fid_log,'\n');
               
        fprintf(fid_log,'Characterization Results:\n');

        fprintf(fid_log,'Read Pairs with >=1 Incorporations:\t%d\n',summary.Incorporation);
        logmatrix(fid_log,summary.IncorporationHistogram,'Incorporation Histogram (0-n)',{'Position' 'Read Pairs'},0:n_cols,'%d');
        
        fprintf(fid_log,'Fraction of Read Pairs with >=1 Incorporations:\t%0.5f\n',summary.Incorporation./numel(P));
        logmatrix(fid_log,summary.IncorporationHistogram_Relative,'Incorporation Histogram, Relative Frequency',{'Position' 'Ratio'},0:n_cols,'%0.5g');
        
        logmatrix(fid_log,summary.TemplateCounts,'Template Base Counts',{'Position' 'A' 'C' 'G' 'U'},1:n_cols,'%d');
        logmatrix(fid_log,summary.TemplateFrequencies,'Template Frequencies',{'Position' 'A' 'C' 'G' 'U'},1:n_cols,'%0.5g');
        logmatrix(fid_log,summary.TemplateFrequenciesVersusControl,'Template Frequencies Versus Control',{'Position' 'A' 'C' 'G' 'U'},1:n_cols,'%0.5g');
        logmatrix(fid_log,summary.TemplateNormalizationFactors,'Template Normalization Factors',{'Position' 'A' 'C' 'G' 'U'},1:n_cols,'%0.5g');
        fprintf(fid_log,'Block Counts Match:\t%d\n',summary.bBlockCountsMatch);
        fprintf(fid_log,'Product Complementary Set N:\t%d\n',summary.ProductComplementarySetN);
        fprintf(fid_log,'Product Mismatch Set N:\t%d\n',summary.ProductMismatchSetN);
        fprintf(fid_log,'Product Unincorporated Set N:\t%d\n',summary.ProductUnincorporatedSetN);
        logmatrix(fid_log,summary.ProductComplementarySetCounts,'Product Complementary Set Base Counts',{'Position' 'A' 'C' 'G' 'U' '-'},1:n_cols,'%d');
        logmatrix(fid_log,summary.ProductUnincorporatedSetCounts,'Product Unincorporated Set Base Counts',{'Position' 'A' 'C' 'G' 'U' '-'},1:n_cols,'%d');
        logmatrix(fid_log,summary.ProductComplementarySetCountsNormalized,'Product Complementary Set Base Counts Normalized',{'Position' 'A' 'C' 'G' 'U'},1:n_cols,'%d');
        logmatrix(fid_log,summary.ProductComplementarySetFrequenciesNormalized,'Product Complementary Set Base Frequencies Normalized',{'Position' 'A' 'C' 'G' 'U'},1:n_cols,'%0.5g');
        logmatrix(fid_log,summary.ProductComplementarySetCountsControlTemplateNormalized,'Product Complementary Set Base Counts Normalized to Control Template',{'Position' 'A' 'C' 'G' 'U'},1:n_cols,'%d');
        logmatrix(fid_log,summary.ProductComplementarySetFrequenciesControlTemplateNormalized,'Product Complementary Set Base Frequencies Normalized to Control Template',{'Position' 'A' 'C' 'G' 'U'},1:n_cols,'%0.5g');
        logmatrix(fid_log,summary.TerminalProductBaseCounts,'Terminal Product Base Counts',{'Position' 'A' 'C' 'G' 'U'},1:n_cols,'%d');
        logmatrix(fid_log,summary.TerminalProductBaseCountsNormalized,'Terminal Product Base Counts Normalized',{'Position' 'A' 'C' 'G' 'U'},1:n_cols,'%d');
        logmatrix(fid_log,summary.TerminalProductBaseFrequenciesNormalized,'Terminal Product Base Frequencies Normalized',{'Position' 'A' 'C' 'G' 'U'},1:n_cols,'%d');
        logmatrix(fid_log,summary.TemplateFirstNullProductBaseCounts,'Template First Null Product Base Counts',{'Position' 'A' 'C' 'G' 'U'},1:n_cols,'%d');
        logmatrix(fid_log,summary.TemplateFirstNullProductBaseCountsNormalized,'Template First Null Product Base Counts Normalized',{'Position' 'A' 'C' 'G' 'U'},1:n_cols,'%d');
        logmatrix(fid_log,summary.TemplateFirstNullProductFrequenciesNormalized,'Template First Null Product Base Frequencies Normalized',{'Position' 'A' 'C' 'G' 'U'},1:n_cols,'%0.5g');
        
        logmatrix(fid_log,summary.TemplateCompSetDimerCounts,'Template Complementary Set Dimer (i, i+1) Counts',{'Dimer' 'AA' 'AC' 'AG' 'AU' 'CA' 'CC' 'CG' 'CU' 'GA' 'GC' 'GG' 'GU' 'UA' 'UC' 'UG' 'UU'},1:(n_cols-1),'%d');
        logmatrix(fid_log,summary.TemplateCompSetDimerFrequencies,'Template Complementary Set Dimer (i, i+1) Frequencies',{'Dimer' 'AA' 'AC' 'AG' 'AU' 'CA' 'CC' 'CG' 'CU' 'GA' 'GC' 'GG' 'GU' 'UA' 'UC' 'UG' 'UU'},1:(n_cols-1),'%0.5g');

        logmatrix(fid_log,summary.TemplateDimerCounts,'Template Dimer (i, i+1) Counts',{'Dimer' 'AA' 'AC' 'AG' 'AU' 'CA' 'CC' 'CG' 'CU' 'GA' 'GC' 'GG' 'GU' 'UA' 'UC' 'UG' 'UU'},1:(n_cols-1),'%d');
        logmatrix(fid_log,summary.TemplateDimerFrequencies,'Template Dimer (i, i+1) Frequencies',{'Dimer' 'AA' 'AC' 'AG' 'AU' 'CA' 'CC' 'CG' 'CU' 'GA' 'GC' 'GG' 'GU' 'UA' 'UC' 'UG' 'UU'},1:(n_cols-1),'%0.5g');

        logmatrix(fid_log,summary.ProductMismatchSetMismatchCounts,'Product Mismatch Set Mismatch Base Counts',{'Position' 'Count'},1:n_cols,'%d');
        logmatrix(fid_log,summary.ProductMismatchSetTerminalMismatchCounts,'Product Mismatch Set Terminal Mismatch Base Counts',{'Position' 'Count'},1:n_cols,'%d');
        logmatrix(fid_log,summary.ProductMismatchSetTerminalMismatchCountsRelativeToAllMismatches,'Product Mismatch Set Terminal Mismatch Counts Relative to All Mismatches',{'Position' 'Count'},1:n_cols,'%0.5g');
        logmatrix(fid_log,summary.ProductMismatchSetMismatchCountsPerProductHistogram,'Product Mismatch Set Mismatch Counts Per Product Histogram',{'N Mismatches per Product' 'Count'},0:n_cols,'%d');
        logmatrix(fid_log,summary.ProductMismatchSetMismatchCountsRelativeToAllIncorporations,'Product Mismatch Set Mismatch Counts Relative to All Incorporations',{'Position' 'Count'},1:n_cols,'%d');
        logmatrix(fid_log,summary.ProductMismatchSetMismatchFollowedByMismatchCounts,'Product Mismatch Set Mismatch Followed by Mismatch Counts',{'Position' 'Count'},1:(n_cols-1),'%d');
        logmatrix(fid_log,summary.ProductMismatchSetMismatchFollowedByNonMismatchCounts,'Product Mismatch Set Mismatch Followed by Non-Mismatch Counts',{'Position' 'Count'},1:(n_cols-1),'%d');
        logmatrix(fid_log,summary.ProductMismatchSetMismatchSequenceSpaceCounts,'Product Mismatch Set Mismatch Sequence Space Counts',{'Position' 'AA' 'AC' 'AG' 'CA' 'CC' 'CU' 'GA' 'GG' 'GU' 'UC' 'UG' 'UU'},1:(n_cols),'%d');
        logmatrix(fid_log,summary.ProductMismatchSetMismatchSequenceSpaceCountsNormalized,'Product Mismatch Set Mismatch Sequence Space Counts Normalized',{'Position' 'AA' 'AC' 'AG' 'CA' 'CC' 'CU' 'GA' 'GG' 'GU' 'UC' 'UG' 'UU'},1:(n_cols),'%d');
        logmatrix(fid_log,summary.ProductMismatchSetMismatchSequenceSpaceFrequenciesNormalized,'Product Mismatch Set Mismatch Sequence Space Frequencies Normalized',{'Position' 'AA' 'AC' 'AG' 'CA' 'CC' 'CU' 'GA' 'GG' 'GU' 'UC' 'UG' 'UU'},1:(n_cols),'%.5g');
        
        % C.7.1 Generate sequence space cube (trimer) plots and data
        seqspace_cube(P,T,options.plots,'logfile',fid_log);
        
        % C.7.2 Generate sequence space cube (trimer) plots and data
        seqspace_cube(P_CompSet,T_CompSet,options.plots,'logfile',fid_log,...
            'product_label','product_CompSet','template_label','template_CompSet');
        
        % C.8 Transition Probabilities. Calculate the position-dependent 
        % base transition probabilities at each position. Generate 
        % transition map plots. At any given position, the next base can 
        % be A, C, G, U, or null (unextended, for product only). Using a 
        % frequentist interpretation of probability, the transition 
        % probability is calculated as the number of transitions observed 
        % normalized by the number of observations at a given position. 
        % This is done in the helper function transition_map(). Counts and 
        % then frequencies are separately computed and visualized for the 
        % product (with and without nulls) and template.
        transition_map(P,T,options.plots,'logfile',fid_log);
        
        % Transition map with complementary set only
        transition_map(P_CompSet,T_CompSet,options.plots,'logfile',fid_log,...
            'description','Transition Map (Complementary Set)',...
            'filenamebase','transition_map_CompSet');
        
        % Close the log file
        fclose(fid_log);
    end    
end
