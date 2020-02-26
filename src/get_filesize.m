function varargout = get_filesize(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = filesize(varargin{:});
else
  [varargout{1}] = filesize(varargin{:});
end
