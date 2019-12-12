function nf = normalization_factor_from_seq(s)
    nf = hash2freq(nt2hash(rna2dna(s)));
    nf = nf(1:4,:);
end
