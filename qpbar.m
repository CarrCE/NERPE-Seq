% QPBAR computes Phred quality score associated with the mean per base 
% error probability given ASCII-encoded quality sequence and an encoding.
%
% Usage: q = qpbar(Q) returns the associated Phred quality score q
% for the quality sequence string Q. By default, it is assumed that the
% quality score is encoded as ASCII characters with value of 
% Phred[i] = Q[i]-33;
%
% Usage: q = qpbar(Q,E) uses the encoding offset E. If E is not specified,
% E=33 is assumed per above.
%
function q_p_bar = qpbar(Q,E)
    if nargin<2, E = 33; end
    % Convert ASCII quality to numeric (phred) quality score
    q_numeric = Q-E;
    % Convert Phred quality score to per base error probability
    p = 10.^(-q_numeric/10);
    % Compute mean per base error probability
    p_bar = mean(p);
    % Compute quality score corresponding to this probability
    q_p_bar = -10*log10(p_bar);
end