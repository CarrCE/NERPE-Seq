function FastqMapSamples(xls)
    [~,~,raw]=xlsread(xls);
    N = size(raw,1)-1;
    for k=1:N
        % Processing sample k of N
        fprintf('FastqMapSamples\n');
        fprintf('Processing sample %d of %d\n',k,N);
        % Get inputs
        fq1 = raw{k+1,1}; fq2 = raw{k+1,2}; index = raw{k+1,3}; out = raw{k+1,4}; 
        template_length = raw{k+1,5}; prefix = raw{k+1,6}; control = raw{k+1,7};
        % Display details
        fprintf('Forward Read File %s\n',fq1);
        fprintf('Reverse Read File %s\n',fq2);
        fprintf('Index Read File   %s\n',index);
        fprintf('Output Directory  %s\n',out);
        fprintf('Template Length   %d\n',template_length);
        fprintf('Prefix            %s\n',prefix);
        fprintf('\n');
        % Make Fastqmaps for input read files
        fastqmap(fq1); fastqmap(fq2); fastqmap(index);
    end
end