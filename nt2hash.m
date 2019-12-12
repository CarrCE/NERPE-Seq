% NT2HASH Converts DNA/RNA nucleotide sequence to a hash code matrix 
% representation M suitable for rapid processing.
%  Usage: M = nt2hash(SEQ)
% SEQ can be a character string sequence or a cell array of sequences.
function [M,hashtable] = nt2hash(seqs)

hashtable = spalloc(255,1,18);
hashtable = NaN*zeros(255,1);
hashtable('-') = 32;    % Gap
hashtable('A') = 1;     % Adenine
hashtable('C') = 2;     % Cytosine
hashtable('G') = 4;     % Guanine
hashtable('T') = 8;     % Thymine
hashtable('U') = 16;    % Uracil
hashtable('R') = sum(hashtable('AG'));  % Purine (A or G)
hashtable('Y') = sum(hashtable('CT'));  % Pyrimidine (C or T)
hashtable('K') = sum(hashtable('GT'));  % Keto (G or T)
hashtable('M') = sum(hashtable('AC'));  % Amino (A or C)
hashtable('S') = sum(hashtable('CG'));  % Strong interaction (3 H bonds)
hashtable('W') = sum(hashtable('AT'));  % Weak interaction (2 H bonds)
hashtable('B') = sum(hashtable('CGT')); % Not A
hashtable('D') = sum(hashtable('AGT')); % Not C
hashtable('H') = sum(hashtable('ACT')); % Not G
hashtable('V') = sum(hashtable('ACG')); % Not T (or U)
hashtable('N') = sum(hashtable('ACGTU')); % Any nucleotide
hashtable('X') = sum(hashtable('ACGTU')); % Any nucleotide

if iscell(seqs),
    numc = numel(seqs);
    M = NaN*ones(numc,size(char(seqs(1)),2));
    for k=1:numc,
        M(k,:) = hashtable(upper(char(seqs(k))));
    end
elseif size(seqs,1)>1,
    M = hashtable(upper(seqs));
else
    M(1,:) = hashtable(upper(seqs));
end

