% 
% FUNCTION OPTIONS = ARGS2OPTIONS(VARARGIN)
%
function options = args2options(varargin)
    tmp = varargin{1};
    if isempty(tmp)
        % no user specified options
        options = [];
    elseif length(tmp)==1
        % assume options struct passed by user
        options = cell2mat(tmp);
    else
        % parse any field-value fields
        odd = mod(1:length(tmp),2);
        bOdd = odd>0;
        f = tmp(bOdd);
        c = tmp(not(bOdd));
        options = cell2struct(c,f,2);
    end
end

