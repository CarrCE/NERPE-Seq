function filter_nerpe(fasta_in,mat_in,fasta_out,mat_out,varargin)
    % filter the fasta file
    filter_nerpe_fasta(fasta_in,fasta_out,varargin{:})
    
    % "filter" the mat file
    load(mat_in)
    
    % Get number of reads 
    [hdr]=read_nerpe_fasta(fasta_out);
    N = numel(hdr);
    
    % add record for how many reads were retained in the fasta file
    summary.FilteredSequencesMatch = N;
    summary.Discard = summary.Discard + (summary.Keep - N);
    summary.Keep = N;
    
    % save the new results
    save(mat_out,'summary','options');
end