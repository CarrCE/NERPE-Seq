% WRITE_NERPE_FASTA writes a fasta-like file of NERPE product-template 
% regions.
%
% write_nerpe_fasta(OUT,H,P,T) writes to the filename OUT the
% product-template regions stored in H, P, T. H is a cell array of headers,
% P is a cell array of product sequences (in 5'-3' orientation) and T is a
% cell array of template sequences (in 3'-5' orientation).
% 
function write_nerpe_fasta(out,H,P,T)
    % Write to file
    fid_fasta = fopen(out,'w');
    N = numel(H);
    for r=1:N
        % Write out fasta as ID with 5' and 3' sequences on
        % separate lines
        fprintf(fid_fasta,'>%s\n',H{r});
        fprintf(fid_fasta,'%s\n',P{r});
        fprintf(fid_fasta,'%s\n',T{r});
    end
    fclose(fid_fasta);
end

            
            