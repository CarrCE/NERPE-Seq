% HASH2FREQ Calculates base pair frequencies for each nucleotide position
% for an array of sequences represented as a hashcode matrix M.
%  Usage: F = hash2freq(M)
%
% F contains frequencies for each base pair position (column) and type of
% nucleotide (row), A, C, G, and T (or U). Use a second output to obtain
% unweighted estimates including for ambiguous nucleotides such as N, 
% meaning A, C, G, or T.
%  Usage: [F,C,U]=hash2freq(M)
%
function [F,C,U] = hash2freq(M)
    Z = zeros(size(M));
    for k=1:32,
        % calculate frequency of kth nucleotide
        % could use sum(bitand(M,k)) but this would
        % multiply count ambiguous nucleotides such as N, 
        % which we will istead take into account as representing
        % 0.25 probability of each nucleotide
        id = find(M==k);
        S = Z;
        S(id) = 1;
        U(k,:) = sum(S,1);
    end
    
    % deal with ambiguous nucleotides
    [tmp,table]=nt2hash('A');
    prob = zeros(32,5);
    % probability  is    [A C G T]
    prob(table('A'),:) = [1 0 0 0 0];
    prob(table('C'),:) = [0 1 0 0 0];
    prob(table('G'),:) = [0 0 1 0 0];
    prob(table('T'),:) = [0 0 0 1 0];
    prob(table('R'),:) = [1 0 1 0 0]*1/2;
    prob(table('Y'),:) = [0 1 0 1 0]*1/2;
    prob(table('K'),:) = [0 0 1 1 0]*1/2;
    prob(table('M'),:) = [1 1 0 0 0]*1/2;
    prob(table('S'),:) = [0 1 1 0 0]*1/2;
    prob(table('W'),:) = [1 0 0 1 0]*1/2;
    prob(table('B'),:) = [0 1 1 1 0]*1/3;
    prob(table('D'),:) = [1 0 1 1 0]*1/3;
    prob(table('H'),:) = [1 1 0 1 0]*1/3;
    prob(table('V'),:) = [1 1 1 0 0]*1/3;
    prob(table('N'),:) = [1 1 1 1 0]*1/4;
    prob(table('X'),:) = [1 1 1 1 0]*1/4;
    prob(table('-'),:) = [0 0 0 0 1];
    
    for k=1:5,
        C(k,:) = prob(:,k)'*U;
    end
    
    for k=1:5,
        F(k,:) = C(k,:)./sum(C,1);
    end
    
end
