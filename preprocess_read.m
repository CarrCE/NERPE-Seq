% PREPROCESS_READ prepares primer extension sequencing data by extracting
% incorporation and template regions for data meeting quality thresholds.
% This helper function works on a single read pair.
%
% preprocess_read(H1,R1,Q1,H2,R2,Q2,H3,R3,Q3) generates results using the 
% default options for forward read with header H1, sequence R1, quality Q1,
% reverse read with header H2, sequence R2, quality Q2, and index read with
% header H3, sequence R3, and quality Q3. The results are returned as a 
% structure.
%
% preprocess_read(...,'property1',value1,'property2',value2,...) generates 
% results using the additional specified options. Valid options are:
%
% Property          Default Value                    Description
% -------------------------------------------------------------------------
%
% PROCESSING
%
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
function results = preprocess_read(H1,R1,Q1,H2,R2,Q2,H3,R3,Q3,varargin)
    %% OPTIONS
    % Set any user specified options
    useroptions = args2options(varargin);
    % Set all options to defaults or user specified options
    options = []; % initial empty options
    options = fieldcheck(options,'q',             30, useroptions);
    options = fieldcheck(options,'qencoding',     33, useroptions);
    options = fieldcheck(options,'fix3','GATCGTCGGACTGTAGAACTCTG',useroptions);
    options = fieldcheck(options,'fix2','AGATCGGAAGAGCACACGTCTGA',useroptions);
    options = fieldcheck(options,'fix1','GCATGCGACTAAACGTCGCATGC', useroptions);
    options = fieldcheck(options,'template_length', 6, useroptions);
    options = fieldcheck(options,'prefix','TCTA', useroptions);
    options = fieldcheck(options,'q_index',       26, useroptions);
    
    % Define output
    results.bHeadersMatch = NaN;
    results.bIndexQualityPass = NaN;
    results.bR1QualityPass = NaN;
    results.bR2QualityPass = NaN;
    results.bQualityPass = NaN;
    results.bR1Fix2Match = NaN;
    results.bR1Fix1Match = NaN;
    results.bR2Fix1Match = NaN;
    results.bR1TemplateLengthMatch = NaN;
    results.bR1PrefixLengthMatch = NaN;
    results.bR2TemplateLengthMatch = NaN;
    results.bR2PrefixLengthMatch = NaN;
    results.bR1PrefixMatch = NaN;
    results.bR1PrefixStart = NaN;
    results.bR2PrefixMatch = NaN;
    results.bR2Fix3Match = NaN;
    results.bProductLengthPass = NaN;
    results.bProductLengthsMatch = NaN;
    results.bSequencesMatch = NaN;
    results.bDiscard = false;
    results.N_incorporations = NaN;
    results.Product = '';
    results.Template = '';
    results.bIncorporation = NaN;
        
    %% B.2.1 Header match
    if ~results.bDiscard
        H1_1 = regexp(H1,' ','split');
        H2_1 = regexp(H2,' ','split');
        results.bHeadersMatch = strcmp(H1_1{1},H2_1{1});
        if ~results.bHeadersMatch
            warning('Headers do not match for read with headers %s and %s',H1,H2);
            results.bDiscard = true;
        end
    end
    % Also demand match with index header H3
    if ~results.bDiscard
        H1_1 = regexp(H1,' ','split');
        H3_1 = regexp(H3,' ','split');
        results.bHeadersMatch = strcmp(H1_1{1},H3_1{1});
        if ~results.bHeadersMatch
            warning('Headers do not match for read with headers %s and %s',H1,H3);
            results.bDiscard = true;
        end
    end
    
    %% B.2.2 Index quality filter; B.2.3 R1 quality filter; B.2.3 R2 quality filter.
    if ~results.bDiscard
        % Check whether Index Read meets quality threshold
        q3_p_bar = qpbar(Q3,options.qencoding);
        results.bIndexQualityPass = q3_p_bar>=options.q_index;
        % Check whether R1 meets quality threshold
        q1_p_bar = qpbar(Q1,options.qencoding);
        results.bR1QualityPass = q1_p_bar>=options.q;
        % Check whether R2 meets quality threshold
        q2_p_bar = qpbar(Q2,options.qencoding);
        results.bR2QualityPass = q2_p_bar>=options.q;
        % Check if we should discard
        results.bQualityPass = and(and(results.bR1QualityPass,results.bR2QualityPass),results.bIndexQualityPass);
        results.bDiscard = ~results.bQualityPass;
    end
    
    %% B.3 Use fix 2 as a quality filter on R1. 
    % Demand that R1 contains a single perfect match to the fix 2 sequence.
    % (R2 does not contain fix 2.)
    if ~results.bDiscard
        R1_fix2_perfect = regexp(R1,options.fix2);
        if or(numel(R1_fix2_perfect)>1,isempty(R1_fix2_perfect))
            % Found either more than one match, or zero matches
            results.bR1Fix2Match = false;
            results.bDiscard = true;
        else
            % Found one match
            results.bR1Fix2Match = true;
        end
    end

    %% B.4 Trim R1. Trim R1 where fix 2 starts 
    % (that is, cut fix 2 and downstream out).
    if ~results.bDiscard
        R1 = R1(1:R1_fix2_perfect-1);
    end
    
    %% B.5 Use fix1 as a quality filter. 
    % Demand that R1 and R2 both contain a perfect match to fix 1.
    if ~results.bDiscard
        % Demand perfect fix1 region in R1
        R1_fix1_perfect = regexp(R1,options.fix1);
        if or(numel(R1_fix1_perfect)>1,isempty(R1_fix1_perfect))
            % Found either more than one match, or zero matches
            results.bR1Fix1Match = false;
            results.bDiscard = true;
        else
            % Found one match
            results.bR1Fix1Match = true;
        end
    end
    if ~results.bDiscard
        % Demand perfect fix2 region in R2
        R2_fix1_perfect = regexp(R2,seqrcomplement(options.fix1));
        if or(numel(R2_fix1_perfect)>1,isempty(R2_fix1_perfect))
            % Found either more than one match, or zero matches
            results.bR2Fix1Match = false;
            results.bDiscard = true;
        else
            % Found one match
            results.bR2Fix1Match = true;
        end    
    end

    %% B.6.1 Confirm R1 Template 
    if ~results.bDiscard
        % Identify template in R1
        R1_template_start = R1_fix1_perfect-options.template_length;
        R1_template_stop = R1_fix1_perfect-1;
        
        % check for incomplete template
        if ~(R1_template_start<1)
            % Complete template
            R1_template = R1(R1_template_start:R1_template_stop);
            results.bR1TemplateLengthMatch=true;
        else
            % Incomplete template
            results.bR1TemplateLengthMatch=false;
            results.bDiscard = true;
        end
    end

    %% B.6.2 Confirm R1 Prefix
    if ~results.bDiscard
        % Identify prefix in R1
        R1_prefix_start = R1_template_start-numel(options.prefix);
        R1_prefix_stop = R1_template_start-1;
        % Check for incomplete prefix
        if ~(R1_prefix_start<1)
            % We have a complete prefix
            results.bR1PrefixLengthMatch = true;
            % Extract complete prefix
            R1_prefix = R1(R1_prefix_start:R1_prefix_stop);
            % Does prefix match?
            if strcmp(R1_prefix,options.prefix)
                results.bR1PrefixMatch=true;
            else
                results.bR1PrefixMatch=false;
                results.bDiscard = true;
            end
        else
            % Incomplete prefix
            results.bR1PrefixLengthMatch = false;
            results.bR1PrefixMatch = false;
            results.bDiscard = true;
        end
    end

    %% B.6.3 Test prefix location in R1 starting at position 1
    if ~results.bDiscard
        % Check prefix starts in R1 at position 1
        if R1_prefix_start == 1
            results.bR1PrefixStart = true;
        else
            results.bR1PrefixStart = false;
            results.bDiscard = true;
        end
    end

    %% B.7.1 Confirm R2 Template 
    if ~results.bDiscard
        % Identify template in R2
        R2_template_start = R2_fix1_perfect+numel(options.fix1);
        R2_template_stop = R2_fix1_perfect+numel(options.fix1)+options.template_length-1;
        % Do we have a complete template?
        if (R2_template_stop<=numel(R2))
            % Template is complete
            results.bR2TemplateLengthMatch = true;
            % Extract the template
            R2_template = seqrcomplement(R2(R2_template_start:R2_template_stop));
        else
            % Template is incomplete
            results.bR2TemplateLengthMatch = false;
            results.bDiscard = true;
        end
    end
    
    %% B.7.2 Confirm R2 Prefix
    if ~results.bDiscard
        % Identify prefix in R2
        R2_prefix_start = R2_template_start+options.template_length;
        R2_prefix_stop = R2_template_start+options.template_length+numel(options.prefix)-1;
        % Do we have a complete prefix?
        if (R2_prefix_stop<=numel(R2))
            % Prefix is complete
            results.bR2PrefixLengthMatch = true;
            % Extract the prefix
            R2_prefix = R2(R2_prefix_start:R2_prefix_stop);
            % Does prefix match?
            if strcmp(seqrcomplement(R2_prefix),options.prefix)
                results.bR2PrefixMatch=true;
            else
                results.bR2PrefixMatch=false;
                results.bDiscard = true;
            end                    
        else
            % Prefix is incomplete
            results.bR2PrefixLengthMatch = false;
            results.bR2PrefixMatch = false;
            results.bDiscard = true;
        end        
    end

    %% B.7.3 Confirm Fix3 match found to right of prefix in R2
    if ~results.bDiscard
        % Identify fix3 in R2
        R2_fix3_start = R2_prefix_stop + 1;
        R2_fix3_stop = R2_fix3_start + numel(options.fix3)-1;
        % Do we have a complete fix3
        if (R2_fix3_stop<=numel(R2))
            % fix3 is complete
            % extract the fix3
            R2_fix3 = R2(R2_fix3_start:R2_fix3_stop);
            % Does fix3 match
            if strcmp(R2_fix3,options.fix3)
                results.bR2Fix3Match=true;
            else
                results.bR2Fix3Match=false;
                results.bDiscard = true;
            end                    
        else
            % fix3 is incomplete
            results.bR2Fix3Match = false;
            results.bDiscard = true;
        end
    end 
    
    %% B.8 Trim R1/R2 in the prefix
    if ~results.bDiscard
        R1 = R1(R1_prefix_stop:end);
        R2 = R2(1:R2_prefix_start);
    end

    %% B.9 Product Identification
    if ~results.bDiscard
        % Determine P1, the product region in R1
        % Note: position is not based on previously determined fix1
        % location because we have now trimmed on the left of R1
        R1_product_start = 1+options.template_length+numel(options.fix1)+1;
        R1_product_stop = numel(R1);
        if R1_product_stop>=R1_product_start
            % Product region has size 1 or greater
            P1 = R1(R1_product_start:R1_product_stop);
        else
            % Product region has 0 length
            P1 = '';
        end
        % Determine P2, the product region in R2
        % Position can be based on previously determined fix1 region
        % because R2 starts with product RC and ends before fix1. We
        % did no left trimming of R2.
        R2_product_start = 1;
        R2_product_stop = R2_fix1_perfect-1;
        if R2_product_stop>=R2_product_start
            % Product region has size 1 or greater
            P2 = seqrcomplement(R2(R2_product_start:R2_product_stop));
        else
            % Product region has 0 length
            P2 = '';
        end
    end

    %% B.10 Product length filter
    if ~results.bDiscard
        maxProductLength = options.template_length+1;
        bR1ProductLengthOk = numel(P1)<=maxProductLength;
        bR2ProductLengthOk = numel(P2)<=maxProductLength;
        if and(bR1ProductLengthOk,bR2ProductLengthOk)
            results.bProductLengthPass = true;
        else
            results.bProductLengthPass = false;
            results.bDiscard = true;
        end
    end
    
    %% B.11 Product equivalent length filter. 
    % Demand that the product length in R1 equal the product length in R2.
    if ~results.bDiscard
        if numel(P1)==numel(P2)
            results.bProductLengthsMatch = true;
        else
            results.bProductLengthsMatch = false;
            results.bDiscard = true;
        end
    end

    %% B.12 Final Sequence Comparison
    if ~results.bDiscard
        if strcmp(R1,seqrcomplement(R2))
            results.bSequencesMatch = true;
        else
            results.bSequencesMatch = false;
            results.bDiscard = true;
        end
    end

    %% B.13 Final processing to yield analyzable data from R1 and conversion to RNA space.
    if ~results.bDiscard
        results.N_incorporations = numel(P1);
        results.bIncorporation = (results.N_incorporations>0);
        results.Product = [dna2rna(P1) repmat('*',[1 options.template_length+1-numel(P1)])];
        results.Template = dna2rna([seqreverse(R1_template) options.prefix(end)]);
    end
end
