% MISMATCH_CONTEXT calculates the context for each mismatch in a set of
% product and template sequences with mismatches. 
%
%
% M = mismatch_context(P,T) returns all mismatches with context
% information.
%
% P is a cell array of product sequences.
% T is a cell array of template sequences.
%
% M is a cell array of size N x 8, where N is the number of
% mismatches in the set of T and P.
%
%
% M = mismatch_context(P,T,b0) allows you to specify the upstream priming
% base on the product strand. By default, the priming base b0='C', the 3'
% end of the nominal fix1 sequence 5''-GCATGCGACTAAACGTCGCATGC-3''
%
%
% M = mismatch_context(P,T,b0,SAVE) writes out table of results to the
% Excel XLSX file specified by SAVE. SAVE should have an '.xlsx' extension.
%
%
% The columns of M are:
% {n score pat0 pat1 pat2 t0 t1 t2}
% 
% n is the position in the product (5'-3') or template (3'-5')
%
% score is an integer from 0 to 5 indicating the pattern of mismatches as
% follows:
%
%   score   pattern                     values of pat0, pat1, pat2
%   0       Correct-Mismatch-Correct    0 m 0
%   1       Correct-Mismatch-Mismatch   0 m m
%   2       Correct-Mismatch-Null       0 m NaN
%   3       Mismatch-Mismatch-Correct   m m 0
%   4       Mismatch-Mismatch-Mismatch  m m m
%   5       Mismatch-Mismatch-Null      m m NaN         
%   
%   m in the above table indicates a mismatch and thus a non-zero value of
%   pat(1), pat(2), or pat(3) as described below. because the table only
%   includes mismatches, the pat1 value is always non-zero.
%
% pat0, pat1, pat2 provides the mismatch identity before the mismatch
% position (0), at the mismatch position (1), and after the mismatch
% position (2). This mismatch identity is encoded as:
%
%   value     0         1  2  3  4  5  6  7  8  9  10 11 12 NaN
%   meaning   correct   AA AC AG CA CC CU GA GG GU UC UG UU X-
% 
%   Meaning is listed as template then product, so pat0 = 3 means the
%   current mismatch has a mismatch at the prior position corresponding to
%   an A in the template and a G in the product. Correctly paired base
%   pairs are assigned a value of 0. Nulls (e.g., no base in the product)
%   are assigned NaN and appear blank in the Excel output.
%
% t0,t1,t2 are the base identities in the template [A, C, G, U] at the 
% template position before the mismatch, at the mismatch, and after 
% the mismatch. Note that these correspond to 5'-3' in product strand, and
% 3'-5' in template strand, so that t0 corresponds to pat0, t1 to pat1, and
% t2 to pat2.
%
% Template identities are encoded within the pat0,pat1,pat2 values except 
% for correct or null bases, where pat0,pat1,pat2 can take on 0 or NaN 
% values. Thus, t0 and t2 are required to reconstruct the precise mismatch 
% context, whereas t1 is redundant but provided for convenience.
% 
% 2020-05-17 Initial version    Christopher E. Carr
% 
function [M] = mismatch_context(P,T,base0,saveas)
    if nargin<4, saveas = ''; end
    if nargin<3
        base0 = 'C';
        warning(['mismatch_context: base at position 0 (base0) ' ...
                 'not specified; assuming priming base C at 3'' ' ...
                 'end of nominal fix1 region '...
                 '5''-GCATGCGACTAAACGTCGCATGC-3''']);
    end

    % determine number of positions N
    N = numel(T{1})-1;
    % number of product and template pairs
    N_MM_set = numel(P);
    mismatch_types = {'AA' 'AC' 'AG' 'CA' 'CC' 'CU' 'GA' 'GG' 'GU' 'UC' 'UG' 'UU'};
    max_mismatches = N_MM_set*N;
    M = cell(max_mismatches,8);
    
    % initialize mismatch counter
    ctr = 0;
    
    % loop through all positions
    for n=1:N
        % Characterize mismatches at this position
        for k=1:N_MM_set
            % Do we have a mismatch at this position
            T_k = T{k}; P_k = P{k};
            T_k_n = T_k(n); P_k_n = P_k(n);
            % If product is not null...
            if P_k_n~='-'
                % Get product-template pair
                pair_k_n = [T_k_n P_k_n];
                % If product-template pair is a mismatch
                bMismatch = strcmpi(mismatch_types,pair_k_n);
                if sum(bMismatch)>0
                    % Get the mismatch index
                    id_mismatch = find(bMismatch);
                    % Extract context for this mismatch
                    if n==1
                        pair_k_nm1 = [seqcomplement(base0) base0];
                    else
                        pair_k_nm1 = [T_k(n-1) P_k(n-1)];
                    end
                    pair_k_np1 = [T_k(n+1) P_k(n+1)];
                    % Is there a mismatch in the n-1 position
                    bMismatch_nm1 = strcmpi(mismatch_types,pair_k_nm1);
                    if sum(bMismatch_nm1)>0
                        id_mismatch_nm1 = find(bMismatch_nm1);
                    else
                        id_mismatch_nm1 = 0;
                    end    
                    % Is there a mismatch in the n+1 position
                    bMismatch_np1 = strcmpi(mismatch_types,pair_k_np1);
                    if sum(bMismatch_np1)>0
                        id_mismatch_np1 = find(bMismatch_np1);
                    else
                        if pair_k_np1(2)=='-'
                            id_mismatch_np1 = NaN;
                        else
                            id_mismatch_np1 = 0;
                        end
                    end    
                    
                    pat = [id_mismatch_nm1 id_mismatch id_mismatch_np1];
                    Tnk = [pair_k_nm1(1) pair_k_n(1) pair_k_np1(1)];
                    
                    % Calculate a score corresponding to a mismatch
                    % pattern. Possible scores are:
                    %
                    % Correct-Mismatch-Correct (0)    0 M 0
                    % Correct-Mismatch-Mismatch (1)   0 M M
                    % Correct-Mismatch-Null (2)       0 M NaN
                    % Mismatch-Mismatch-Correct (3)   M M 0
                    % Mismatch-Mismatch-Mismatch (4)  M M M
                    % Mismatch-Mismatch-Null (5)      M M NaN
                    score = (pat(1)>0)*3 + max((pat(3)>0),isnan(pat(3))*2);
                    
                    % Pattern is
                    % 0         1    2    3    4    5    6   7     8   9    10   11   12
                    % correct {'AA' 'AC' 'AG' 'CA' 'CC' 'CU' 'GA' 'GG' 'GU' 'UC' 'UG' 'UU'};
                    % mismatches are 'TP', e.g., template then product base
                    
                    % Store results
                    ctr = ctr + 1;
                    M{ctr,1} = n;
                    M{ctr,2} = score;
                    M{ctr,3} = pat(1);
                    M{ctr,4} = pat(2);
                    M{ctr,5} = pat(3);
                    M{ctr,6} = Tnk(1);
                    M{ctr,7} = Tnk(2);
                    M{ctr,8} = Tnk(3);
                    
                    % Mutation context is represented as
                    % {n score pat(1) pat(2) pat(3) Tnk(1) Tnk(2) Tnk(3)}
                end
            end
        end
    end
    % Remove any extra entries in the mutation context matrix by saving 
    % only the first ctr rows of the matrix
    M = M(1:ctr,:);
    
    % Write out file
    if ~isempty(saveas)
        M_table = cell2table(M,'VariableNames',{'position','score','pattern0','pattern1','pattern2','Template0','Template1','Template2'});
        writetable(M_table,saveas);
    end
end
