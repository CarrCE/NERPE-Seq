%
% FUNCTION S = FIELDCHECK(S,FIELD,DEFAULT)
%
function s = fieldcheck(s,field,default,useroptions)
    % allow user value, if present
    % to override default option
    if isfield(useroptions,field),
        default = useroptions.(field);
    end
    % fill in default option
    if isempty(s),
        s= struct(field,default);
    elseif ~isfield(s,field),
        s.(field) = default;
    end
end