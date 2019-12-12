% FASTQMAP Produces a fastq file index useful in fast reading of a large fastq file
%  Usage: fastamap(fastafile)
%  The index file will be saved to [fastafile '.map.mat']

%
% (C) Christopher E. Carr (chrisc@mit.edu) 2011/06/02
%
% USE AT YOUR OWN RISK. THIS ITEM IS PROVIDED "AS IS" WITHOUT WARRANTY OF
% ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO
% THE QUALITY AND PERFORMANCE OF THE ITEM IS WITH YOU. SHOULD THE PROGRAM
% PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL ERRORS.
%
% IN NO EVENT WILL ANY COPYRIGHT HOLDER OR ANY OTHER PARTY WHO MAY MODIFY
% AND/OR REDISTRIBUTE THE PROGRAM, BE LIABLE TO YOU FOR DAMAGES, INCLUDING
% ANY GENERAL, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF
% THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO
% LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU
% OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER
% PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGES.
%

function [fn,map,hdr] = fastqmap(fastafile,N,overwrite)
    if nargin<2, N=1; end
    if nargin<3, overwrite=0; end

    % Check for input file
    if ~exist(fastafile,'file'), error('File %s not found.',fastafile); end
    % Check for existence of output file
    [p,f,x]=fileparts(fastafile);
    fn = fullfile(p,[f x '.map.mat']);
    if or(~exist(fn,'file'),overwrite)
        % Make map
        % preallocate memory
        map = NaN(N,1);
        hdr = repmat({''},[N 1]);

        % open the file
        fid = fopen(fastafile,'r');

        bDone = 0;
        id = 0;
        maxpos = NaN; % useful for testing purposes to limit execution time

        while (~bDone),
            pos = ftell(fid);
            line = fgetl(fid);
            if isnumeric(line),
                % found end of file
                bDone = 1;
            elseif isempty(line),
                % do nothing
            else
                if strcmp(line(1),'@'),
                    % found header
                    id = id + 1;
                    map(id,1) = pos;
                    hdr(id,1) = {line(2:(numel(line)))};
                    line = fgetl(fid); % sequence
                    line = fgetl(fid); % quality header
                    line = fgetl(fid); % quality data
                    %chk(id,1) = uint16(sum(line(2:(numel(line)-2))));
                    % Typical header:
                    % @HWI-ST728:7:1101:1438:2019#0/1
                    % chk is checksum on "HWI-ST728:7:1101:1438:2019#0"

                    % debug: show header and seq
                    % disp(line);
                    % pause;
                    if ~mod(id,250000), fprintf('%s processed %d reads\n',fn,id); end
                end
            end

            if ~isnan(maxpos)
                if pos>=maxpos
                    bDone = 1;
                end
            end
        end

        save(fn,'map');
        fclose(fid);
        fprintf('fastqmap made for %s; processed %d reads.\n',fn,id);
    else
        load(fn); N = numel(map);
        fprintf('fastqmap already exists for %s; found %d reads.\n',fn,N);
    end
end