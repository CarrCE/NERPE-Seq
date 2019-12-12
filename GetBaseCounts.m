function C = GetBaseCounts(Seqs,maxbases,map,maxlen)
    C = zeros(maxbases,maxlen); col_id = 1:maxlen;
    for k=1:numel(Seqs)
        n_k = numel(Seqs{k});
        idx_k = sub2ind([maxbases maxlen],full(map(Seqs{k})),col_id(1:n_k));
        C(idx_k) = C(idx_k)+1;
    end
end
