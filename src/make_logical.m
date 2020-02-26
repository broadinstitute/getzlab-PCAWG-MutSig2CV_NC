function S = make_logical(S,varargin)
%
% for a given struct S,
%   performs str2double on all the fields specified and then converts them to logical (boolean)
%
% Mike Lawrence 2013-02-13

fields = {};
for i=1:length(varargin)
  fields = [fields varargin{i}];
end

for i=1:length(fields)
  if isfield(S,fields{i})    
    x = getfield(S,fields{i});
    if ~islogical(x) && ~isnumeric(x)
      x = str2double(x);
    end
    if ~islogical(x)
      x = (x~=0);
    end
    S = setfield(S,fields{i},x);
  else
    fprintf('No such field: %s\n', fields{i});
  end
end
