% TRANSITION_MAP calculates and visualizes transition probabilities for
% NERPE blocks, e.g. matched product and template regions.
%
% This is a helper function called from CHARACTERIZE.m
% See also C.8 Transition Probabilities.
%
function transition_map(P,T,outfolder,varargin)    
    %% OPTIONS
    % Set any user specified options
    useroptions = args2options(varargin);
    % Set all options to defaults or user specified options
    options = []; % initial empty options
    options = fieldcheck(options,'fontsize', 12, useroptions);
    options = fieldcheck(options,'figformat', {'eps' '-depsc'}, useroptions);
    options = fieldcheck(options,'fontname', 'Helvetica', useroptions);
    options = fieldcheck(options,'filenamebase','transition_map',useroptions);
    options = fieldcheck(options,'description','Transition Map',useroptions);
    options = fieldcheck(options,'logfile', '', useroptions);
    
    % make output folder if necessary
    if ~exist(outfolder,'dir'), mkdir(outfolder); end

    % Get template length
    t_len = numel(T{1})-1;
    
    % Write out results
    if ~isempty(options.logfile)
        % Logfile is either a filename or a file identifier
        % Check for file identifier 
        if ischar(options.logfile)
            fid = fopen(filename,'w');
        elseif (fopen(options.logfile)~=-1)
            % Have valid file handle
            fid = options.logfile;
        else
            fid = NaN;
        end
    else
        fid = NaN;
    end

    % Get base transition frequencies for products
    s = '-ACGU';
    C_P = transition_counts(P,s);
    % Calculate frequency
    F_P = C_P./repmat(sum(C_P,1),numel(s),1);

    % Get base transition frequencies for template
    s = 'ACGU';
    C_T = transition_counts(T,s);
    % Calculate frequency
    F_T = C_T./repmat(sum(C_T,1),numel(s),1);

    % Get colormap
    cmap = parula;

    % Generate a heatmap, NOT USING MATLAB BUILT IN FUNCTION (BAD)
    % https://www.mathworks.com/help/matlab/ref/heatmap.html
    % Heatmap does not support custom tick labels.

    % Labels for template
    ytl = {'A' 'C' 'G' 'U'};
    xtl = repmat(ytl,[1 t_len]);

    % Reoriented version
    figure('color',[1 1 1]);
    set(gcf, 'Renderer', 'Painters');
    imagesc(F_T');
    colormap(cmap);
    set(gca,'ytick',1:24,'yticklabel',xtl);
    set(gca,'xtick',1:4,'xticklabel',ytl(1:4));
    colorbar;
    ylabel('Current Base');
    xlabel('Next Base');
    title('Template Transition Frequencies');

    % Save figure
    fn = [options.filenamebase '_template'];
    print(fullfile(outfolder,[fn '.' options.figformat{1}]),options.figformat{2});
    savefig(fullfile(outfolder,[fn '.fig']));
    close(gcf);

    % Labels for product
    ytl = {'-' 'A' 'C' 'G' 'U'};
    xtl = repmat(ytl,[1 t_len]);

    % Reoriented version
    % With nulls
    figure('color',[1 1 1]);
    set(gcf, 'Renderer', 'Painters');
    imagesc(F_P');
    colormap(cmap);
    set(gca,'ytick',1:30,'yticklabel',xtl);
    set(gca,'xtick',1:5,'xticklabel',ytl(1:5));
    colorbar;
    ylabel('Current Base');
    xlabel('Next Base');
    title('Product Transition Frequencies');

    % Save figure
    fn = [options.filenamebase '_product'];
    print(fullfile(outfolder,[fn '.' options.figformat{1}]),options.figformat{2});
    savefig(fullfile(outfolder,[fn '.fig']));
    close(gcf);

    % Reoriented version
    % Without nulls
    idx=(find(~strcmpi(xtl,'-')));
    figure('color',[1 1 1]);
    set(gcf, 'Renderer', 'Painters');
    imagesc(F_P(2:end,idx)');
    colormap(cmap);
    set(gca,'ytick',1:sum(idx),'yticklabel',xtl(idx));
    set(gca,'xtick',1:4,'xticklabel',ytl(2:5));
    colorbar;
    ylabel('Current Base');
    xlabel('Next Base');
    title('Product Transition Frequencies (no nulls)');

    % Save figure    
    fn = [options.filenamebase '_product_no_nulls'];
    print(fullfile(outfolder,[fn '.' options.figformat{1}]),options.figformat{2});
    savefig(fullfile(outfolder,[fn '.fig']));
    close(gcf);
    
    % If appropriate, log to file
    if ~isnan(fid)
        % Have valid file handle
        % Print data to the file
        fprintf(fid,sprintf('%s:\n',options.description));
        logmatrix(fid,C_P','Product Transition Counts',{'Current Base (Row)/Next Base (Col)' '1-','1A','1C','1G','1U','2-','2A','2C','2G','2U','3-','3A','3C','3G','3U','4-','4A','4C','4G','4U','5-','5A','5C','5G','5U','6-','6A','6C','6G','6U'},{'-' 'A' 'C' 'G' 'U'},'%d');
        fprintf(fid,'\n');
        logmatrix(fid,F_P','Product Transition Frequencies',{'Current Base (Row)/Next Base (Col)','1-','1A','1C','1G','1U','2-','2A','2C','2G','2U','3-','3A','3C','3G','3U','4-','4A','4C','4G','4U','5-','5A','5C','5G','5U','6-','6A','6C','6G','6U'},{'-' 'A' 'C' 'G' 'U'},'%0.4f');
        fprintf(fid,'\n');
        logmatrix(fid,C_T','Template Transition Counts',{'Current Base (Row)/Next Base (Col)','1A','1C','1G','1U','2A','2C','2G','2U','3A','3C','3G','3U','4A','4C','4G','4U','5A','5C','5G','5U','6A','6C','6G','6U'},{'A' 'C' 'G' 'U'},'%d');
        fprintf(fid,'\n');
        logmatrix(fid,F_T','Template Transition Frequencies',{'Current Base (Row)/Next Base (Col)','1A','1C','1G','1U','2A','2C','2G','2U','3A','3C','3G','3U','4A','4C','4G','4U','5A','5C','5G','5U','6A','6C','6G','6U'},{'A' 'C' 'G' 'U'},'%0.4f');
        fprintf(fid,'\n');
    end
end
