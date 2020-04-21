% filter_nepe_fasta(in,out,'defined','GCCNNN')
% filter_nepe_fasta(in,out,'gc',0.80)

function filter_nerpe_fasta(in,out,varargin)
    %% OPTIONS
    % Set any user specified options
    useroptions = args2options(varargin);
    % Set all options to defaults or user specified options
    options = []; % initial empty options
    options = fieldcheck(options,'defined','', useroptions);
    options = fieldcheck(options,'gc', NaN, useroptions);

    if ~isnan(options.gc)
        % Run in GC mode
        % Require input file to exist
        assert(exist(in,'file')==2);
        % Open the file
        [H,P,T]=read_nerpe_fasta(in,4E5);
        % Filter the file
        N = numel(H); b = false(N,1);
        parfor k=1:N
            bc_k = basecount(T{k});
            gc_k = (bc_k.C + bc_k.G)/(bc_k.A + bc_k.C+bc_k.G+bc_k.T);
            if gc_k>=options.gc 
                b(k)=1;
            end
        end
        % Write the output file
        write_nepe_fasta(out,H(b),P(b),T(b));
        
    elseif ~isempty(options.defined)
        % Run in defined mode
        % Require input file to exist
        assert(exist(in,'file')==2);
        % Open the file
        [H,P,T]=read_nerpe_fasta(in,4E5);
        % Convert query sequence to regular expression
        q_regexp = strrep(seq2regexp(options.defined,'Ambiguous',false),'T','U');
        % Filter the file
        N = numel(H); b = false(N,1);
        for k=1:N
            % Version 1.0
            % idx = regexp(rna2dna(T{k}),q_regexp);
            % Version 1.1 modified to filter in RNA space
            idx = regexp(T{k},q_regexp);
            if ~isempty(idx)
                % Found defined sequence at position 1
                b(k)=(idx==1);
            end
        end
        % Write the output file
        write_nerpe_fasta(out,H(b),P(b),T(b));
    else
        % Nothing specified
        warning('filter_nerpe_fasta: no valid processing option specified');
    end
end

