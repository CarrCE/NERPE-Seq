% TRANSITION_COUNTS calculates the symbol transition counts for
% a set of sequences, from which transition frequencies or probabilities
% can easily be estimated.
%
% C = transition_counts(SEQS) returns C, a k x N matrix of
% counts, where each row 1...k corresponds to a particular symbol and
% each column 1...N corresponds to a sequence position. SEQS is a
% cell array of sequences. All sequences are assumed to be the same length.
% 
% [C,S]=transition_counts(SEQS) also returns the set of symbols S, a
% column vector for which each row corresponds to a unique symbol in SEQS
% and the corresponding row in C.
% 
% C = transition_counts(SEQS,S) forces use of the symbols S, useful to also
% count zeros when certain sequences may not be represented in the
% sequences SEQS.
%
function [C,s] = transition_counts(seqs,s)
    N = numel(seqs{1});
    if nargin<2, s = unique([seqs{:}]); end
    N_s = numel(s);
    % Symbol map
    map = full(sparse(ones(size(s)),double(s),1:N_s));
    % Preallocate counts matrix
    C = zeros(numel(s),numel(s)*(N-1));
    % Look through all positions
    for k=1:(N-1)
        % Get dimers corresponding to position k and position k+1
        dimers = cellfun(@(x)([x(k) x(k+1)]),seqs,'UniformOutput',false);
        offset = N_s * (k-1);
        tmp = cellfun(@(x)(map(x)'),dimers,'UniformOutput',false);
        cr = [tmp{:}]';
        [u,ia,ic] = unique(cr,'rows');
        counts_k = accumarray(ic,1);
        % Put these counts into C
        C(sub2ind(size(C),u(:,2),u(:,1)+offset))=counts_k;
    end
end