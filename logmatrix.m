function logmatrix(fid,M,title,rows,cols,format)
    % Title
    if ~isempty(title)
        fprintf(fid,'%s\n');
    end
    % Title
    fprintf(fid,'%s:\n',title);
    % Cols
    fprintf(fid,'%s\t',rows{1}); % For row headers
    for k=1:numel(cols)
        if isnumeric(cols)
            fprintf(fid,'%d',cols(k));
        elseif iscell(cols)
            fprintf(fid,'%s',cols{k});
        else
            % not handled
        end
        if k<numel(cols)
            fprintf(fid,'\t');
        else
            fprintf(fid,'\n'); 
        end
    end
    % Rows
    for k=2:numel(rows)
        % Row heading
        fprintf(fid,'%s\t',rows{k});
        % Row data
        row_tmp = sprintf([format '\t'],M(k-1,:));
        fprintf(fid,'%s\n',row_tmp(1:end-1));    
    end
    % Blank row for spacing
    fprintf(fid,'\n');
end
