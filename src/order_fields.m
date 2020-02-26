function varargout = order_fields(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = orderfields(varargin{:});
else
  [varargout{1}] = orderfields(varargin{:});
end
