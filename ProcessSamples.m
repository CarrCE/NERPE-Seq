% PROCESSSAMPLES performed NERPE-Seq analysis on a set of samples.
%
% ProcessSamples(XLS) performs NERPE-Seq analysis on the set of samples
% described by the rows of the Microsoft Excel file associated with
% filename XLS. This file should contain a header row with the headings as
% follows below, with one row per sample:
%
% Column 1: Heading 'fq1'
% Description: forward read filename
% Example row entry: '../../data/Run_F/6NC4_S3_L001_R1_001.fastq'
%
% Column 2: Heading 'fq2'
% Description: reverse read filename
% Example row entry: '../../data/Run_F/6NC4_S3_L001_R2_001.fastq'
%
% Column 3: Heading 'index'
% Description: index read filename
% Example row entry: '../../data/Run_F/6NC4_S3_L001_I1_001.fastq'
%
% Column 4: Heading 'out'
% Description: output folder
% Example row entry: './out/6NC4'
%
% Column 5: Heading 'template_length'
% Description: length of template sequence
% Example row entry: 6
%
% Column 6: Heading 'prefix'
% Description: sequence in DNA space {A, C, G, T} of prefix
% Example row entry: 'TCTA'
%
% Column 7: Heading 'control'
% Description: MAT file of any associated control condition
% Example row entry: '' or './out/Control/Characterized.mat'
%
% Column 8: Heading 'normalization_factor'
% Description: The normalization factor is a matrix of base frequencies
% with four rows, corresponding to A, C, G, U, and one column for each
% template position. This can be specified as a matrix of values
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
function ProcessSamples(xls)
    [~,~,raw]=xlsread(xls);
    N = size(raw,1)-1;
    for k=1:N
        % Processing sample k of N
        fprintf('Processing sample %d of %d\n',k,N);
        % Get inputs
        fq1 = raw{k+1,1}; fq2 = raw{k+1,2}; index = raw{k+1,3}; out = raw{k+1,4}; 
        template_length = raw{k+1,5}; prefix = raw{k+1,6}; control = raw{k+1,7};
        normalization_factor = raw{k+1,8};
        % Display details
        fprintf('Forward Read File %s\n',fq1);
        fprintf('Reverse Read File %s\n',fq2);
        fprintf('Index Read File   %s\n',index);
        fprintf('Output Directory  %s\n',out);
        fprintf('Template Length   %d\n',template_length);
        fprintf('Prefix            %s\n',prefix);
        fprintf('\n');
        % Make Fastqmaps for input read files
        fastqmap(fq1); fastqmap(fq2); fastqmap(index);
        % Do PreProcessing
        preprocess(fq1,fq2,index,out,'template_length',template_length,'prefix',prefix);
        % Do Characterization
        % check for normalization_factor as expression or as string
        if ~contains(normalization_factor,'[')
            nf = normalization_factor_from_seq(normalization_factor);
        else
            nf = eval(normalization_factor);
        end
        characterize(out,'control',control,'normalization_factor',nf);
    end
end