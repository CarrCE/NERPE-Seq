%READ_NERPE_FASTA reads in a NERPE-style fasta file with entries similar to
% >HEADER
% PRODUCT
% TEMPLATE
%
% [H,P,T]=read_nerpe_fasta(fasta), reads contents of fasta file with NEPE
% formatting and returns a cell array of headers H, products P, and
% templates T.
%
function [H,P,T]=read_nerpe_fasta(fasta,preallocate)
    % Preallocate memory
    if nargin<2, preallocate = 1; end
    H = cell(preallocate,1);
    P = cell(preallocate,1);
    T = cell(preallocate,1);
    % Open file
    fid = fopen(fasta,'r');
    % Prepare to read in data
    bDone = false; ctr = 0;
    while ~bDone
        % Read a record
        l1 = fgetl(fid);
        l2 = fgetl(fid);
        l3 = fgetl(fid);
        % Check for EOF
        if ~ischar(l1)
            % EOF (fgetl returns -1)
            bDone = true;
        else
            % Verify record
            if (l1(1) == ">")
                ctr = ctr + 1;
                H(ctr) = {l1(2:end)};
                P(ctr) = {l2};
                T(ctr) = {l3};
            else
                % Error
                fprintf('Error');
            end
        end
    end
    % Restrict output to actual reads found 
    % (in case we preallocated too much memory)
    H = H(1:ctr);
    P = P(1:ctr);
    T = T(1:ctr);
end