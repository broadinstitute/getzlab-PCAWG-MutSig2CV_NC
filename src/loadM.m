function M = loadM(indir,flds)

if (nargin~=1 && nargin~=2) || nargout~=1 || ~ischar(indir), error('usage:  M = loadM(indir)'); end

if ~exist(indir,'dir')
  error('directory %s does not exist',indir);
end

M = [];
d=[];
d.file = direc([indir '/*']);
d.is_mat = grepm('\.mat$',d.file);
d.datenum = get_filedatenum(d.file); d = sort_struct(d,'datenum');
d = parse_in(d,'file',[indir '/([^\.]*)(\.mat)?$'],'field');
for i=1:slength(d)
  if exist('flds','var') && ~ismember(d.field{i},flds), continue; end
  clear tmp;
  if d.is_mat(i)
    if exist(d.file{i},'file')
      load(d.file{i},'tmp');
      if ~exist('tmp','var'), error('%s does not contain "tmp" variable',d.file{i}); end
    else
      error('error loading %s',d.file{i});
    end
  elseif exist(d.file{i},'dir')
    tmp = load2(d.file{i});
  else
    error('members should be structdir or matfile');
  end
  M = setfield(M,d.field{i},tmp);
end
