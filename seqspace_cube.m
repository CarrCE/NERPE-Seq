function seqspace_cube(P,T,outfolder,varargin)
% Extract trimers at a specified position

% When analyzing data, you would typically compare product sequence space
% to the template sequence space of a control. In other words, run this
% function for both a specific sample and the relevant control.

    %% OPTIONS
    % Set any user specified options
    useroptions = args2options(varargin);
    % Set all options to defaults or user specified options
    options = []; % initial empty options
    options = fieldcheck(options,'fontsize', 12, useroptions);
    options = fieldcheck(options,'figformat', {'eps' '-depsc' '-painters'}, useroptions);
    options = fieldcheck(options,'fontname', 'Helvetica', useroptions);
    options = fieldcheck(options,'logfile', '', useroptions);
    options = fieldcheck(options,'product_label','product',useroptions);
    options = fieldcheck(options,'template_label','template',useroptions);
    options = fieldcheck(options,'guideline',true,useroptions);
    
    % make output folder if necessary
    if ~exist(outfolder,'dir'), mkdir(outfolder); end

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

    % Define way to map RNA ASCII to RNA base
    bases = 'ACGU';
    ascii_2_baseid = full(sparse([1 1 1 1],double(bases),[1 2 3 4]));
    
    % Define triplet combos
    combos = permn(1:4,3);
    % Define triplet combos of bases
    % base_combos = bases(combos);
    
    % Remove any * or - bases in product
    P = strrep(P,'-','');
    P = strrep(P,'*','');
    p_len = cellfun(@numel,P);

    % Get template length
    t_len = numel(T{1})-1;
    
    % Extract trimers at a specified position
    for k=1:(t_len-3)
        % Select Product Length for Sequence Space Analysis
        bMatch = (p_len>=k+2);
        P_Match = P(bMatch);
        P_Match3 = cellfun(@(x)(x(k:k+2)),P_Match,'UniformOutput',false);        

        % Map to coordinates
        tmp = cellfun(@(x)(ascii_2_baseid(x)),P_Match3,'UniformOutput',false);
        xyz = reshape([tmp{:}],3,numel(tmp));
        
        figure('color',[1 1 1]);        
        [u,~,ic]=unique(xyz','rows','stable');
        counts = accumarray(ic, 1); % Count Occurrences
        c = counts/sum(counts)*10000;
        for j=1:numel(c)
            plot3(u(j,1),u(j,2),u(j,3),'ok','MarkerSize',2*(c(j)+1)^(1/3),'MarkerFaceColor','k'); hold on; % Unfilled
        end        
                
        grid on;
        xlabel('Base 1'); ylabel('Base 2'); zlabel('Base 3');
        set(gca,'xtick',[1 2 3 4],'xticklabel',{'A' 'C' 'G' 'U'});
        set(gca,'ytick',[1 2 3 4],'yticklabel',{'A' 'C' 'G' 'U'});
        set(gca,'ztick',[1 2 3 4],'zticklabel',{'A' 'C' 'G' 'U'});
        set(gca,'xlim',[1 4],'ylim',[1 4],'zlim',[1 4]);
        set(gca,'fontname',options.fontname,'fontsize',options.fontsize);
        
        % add guide line
        if options.guideline
            for jj=1:4
                for kk=1:4
                    hl = plot3([jj jj],[1 4],[kk kk],'-');
                    set(hl,'Color',[0 0 0 0.3]);
                end
            end
        end
                
        camproj('perspective');
        set(gca,'CameraPosition',[-1.61 -22.6 7.9]);
        
        fn = sprintf('seqspace_cube_%s_k%d',options.product_label,k);
        print(options.figformat{3},fullfile(outfolder,[fn '.' options.figformat{1}]),options.figformat{2});
        savefig(fullfile(outfolder,[fn '.fig']));
        
        % If appropriate, log to file
        if ~isnan(fid)
            % Have valid file handle
            % Print data to the file
            fprintf(fid,'Sequence Space Cube (Trimer Frequencies): Product, Position %d\n',k);
            % For each combo, check for count
            for j=1:size(combos,1)
                % Check if this combo exists in variable count
                id = find(ismember(u,combos(j,:),'rows'));
                if ~isempty(id)
                    % Count exists, write out base_combo and count
                    fprintf(fid,'%s\t%d\n',bases(combos(j,:)),counts(id));
                else
                    % Count doesn't exist, write out base_combo and zero
                    fprintf(fid,'%s\t%d\n',bases(combos(j,:)),0);
                end
            end
            fprintf(fid,'\n');
        end
        
    end
 
    % Should be no gaps in the template
    t_len = cellfun(@numel,T);

    for k=1:(t_len-3)
        % Select Template Length for Sequence Space Analysis
        bMatch = (t_len>=k+2);
        T_Match = T(bMatch);
        T_Match3 = cellfun(@(x)(x(k:k+2)),T_Match,'UniformOutput',false);        

        % Map to coordinates
        tmp = cellfun(@(x)(ascii_2_baseid(x)),T_Match3,'UniformOutput',false);
        xyz = reshape([tmp{:}],3,numel(tmp));

        figure('color',[1 1 1]);
        [u,~,ic]=unique(xyz','rows','stable');
        counts = accumarray(ic, 1); % Count Occurrences
        c = counts/sum(counts)*10000;
        for j=1:numel(c)
            plot3(u(j,1),u(j,2),u(j,3),'ok','MarkerSize',2*(c(j)+1)^(1/3),'MarkerFaceColor','k'); hold on; % Unfilled
        end
        
        grid on;
        xlabel('Base 1'); ylabel('Base 2'); zlabel('Base 3');
        set(gca,'xtick',[1 2 3 4],'xticklabel',{'A' 'C' 'G' 'U'});
        set(gca,'ytick',[1 2 3 4],'yticklabel',{'A' 'C' 'G' 'U'});
        set(gca,'ztick',[1 2 3 4],'zticklabel',{'A' 'C' 'G' 'U'});
        set(gca,'xlim',[1 4],'ylim',[1 4],'zlim',[1 4]);
        set(gca,'fontname',options.fontname,'fontsize',options.fontsize);
        
        % add guide line
        if options.guideline
            for jj=1:4
                for kk=1:4
                    hl = plot3([jj jj],[1 4],[kk kk],'-');
                    set(hl,'Color',[0 0 0 0.3]);
                end
            end
        end
        
        camproj('perspective');
        set(gca,'CameraPosition',[-1.61 -22.6 7.9]);

        fn = sprintf('seqspace_cube_%s_k%d',options.template_label,k);
        print(options.figformat{3},fullfile(outfolder,[fn '.' options.figformat{1}]),options.figformat{2});
        savefig(fullfile(outfolder,[fn '.fig']));
        
        % If appropriate, log to file
        if ~isnan(fid)
            % Have valid file handle
            % Print data to the file
            fprintf(fid,'Sequence Space Cube (Trimer Frequencies): Template, Position %d\n',k);
            % For each combo, check for count
            for j=1:size(combos,1)
                % Check if this combo exists in variable count
                id = find(ismember(u,combos(j,:),'rows'));
                if ~isempty(id)
                    % Count exists, write out base_combo and count
                    fprintf(fid,'%s\t%d\n',bases(combos(j,:)),counts(id));
                else
                    % Count doesn't exist, write out base_combo and zero
                    fprintf(fid,'%s\t%d\n',bases(combos(j,:)),0);
                end
            end
            fprintf(fid,'\n');
        end
        
    end
end
