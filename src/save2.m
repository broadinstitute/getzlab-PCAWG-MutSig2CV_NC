function save2(X,outdir)

if nargin~=2 || ischar(X) || ~ischar(outdir), error('usage:  save2(X,outdir)'); end
if ~isstruct(X), error('X should be a struct'); end

if exist(outdir,'file') && ~exist(outdir,'dir')
  error('%s already exists and is not a directory',outdir);
end

ede(outdir);

f = fieldnames(X);
nf = length(f);

%diff fields in struct/save directory
D = [];
D.dirs = direc(outdir);
D = reorder_struct(D, ~grepm('md5$', D.dirs));
D.fields = regexprep(D.dirs, '.*/(.*)\.(mat|save2_by_chunks)', '$1');
D.stem = regexprep(D.dirs, '(.*)\.(mat|save2_by_chunks|md5)', '$1');
D.is_dir = grepm('save2_by_chunks$', D.dirs);

for i = find(~ismember(D.fields, f))',
  fprintf('Removing extraneous saved field "%s" ...\n', D.fields{i});

  if D.is_dir(i), rmdir(D.dirs{i}); end
  delete([D.stem{i} '.*'])
end

for i=1:nf,if nf>1, fprintf('%d/%d ',i,nf); end
  t1 = tic;
  tmp = getfield(X,f{i});

  if isnumeric(tmp) && ~isreal(tmp)
    fprintf('WARNING: numerical field "%s" is marked as being complex.\n');
    fprintf('  This is usally due to a bug in str2doubleq().  Use str2doubleq_wrapper().\n');
    fprintf('  Converting to real so that md5 checksum can be computed.\n');
    tmp = real(tmp);
  end
 
  %calculate checksum for this field; if it's unperturbed, then skip it

  if iscell(tmp),
    cs = md5(tmp, 'Array');
  else
    cs = md5(tmp);
  end

  md5path = [outdir '/' f{i} '.md5'];
  
  if exist(md5path) == 2,
    cs_md5 = load_lines(md5path);
    if length(cs_md5) ~= 1, error('Error loading checksum!'); end
    cs_md5 = decell(cs_md5);
    if length(cs_md5) ~= 32, error('Malformed checksum!'); end
    
    if strcmp(cs_md5, cs), 
      fprintf('[skip] ', f{i});
      continue;
    else
      fprintf('(%s) ', f{i});
    end
  end

  % check to see if it's too large to save with "save"
  w = whos('tmp');
  sz = w.bytes;
  save_in_chunks = (sz>=2e9);

  if ~save_in_chunks
    fname = [outdir '/' f{i} '.mat'];
    save(fname,'tmp');
    filesz = get_filesize(fname);
    if filesz==128 % it didn't work
      delete(fname);
      save_in_chunks = true;
    end
  end

  if save_in_chunks
    try
      tmp2 = tmp(1:2);
    catch me
      fprintf('Not saving: %s\n',f{i});
      save_in_chunks = false;
    end
  end

  if save_in_chunks
    n = slength(X);
    nchunks = ceil(sz/1e8);
    chunksize = ceil(n/nchunks);
    outdir2 = [outdir '/' f{i} '.save2_by_chunks'];
    ede(outdir2);
    ci=0;
    for i=1:chunksize:n
      ci=ci+1;
      tmp2 = tmp(i:min(n,i+chunksize-1), :, :, :, :, :, :, :, :, :, :, :, :, :, :, :);
      fname2 = sprintf([outdir2 '/chunk%08d'],ci);
      save(fname2,'tmp2');
      filesz = get_filesize(fname2);
      if filesz==128 % it didn't work
        error('failed to save in chunks');
      end
    end
  end

  %save checksum
  save_lines(cellstr(cs), md5path);

  t2 = toc(t1);
  if t2<1.05, pause(1.05-t2); end   % so that fields can be reload2d in the save order they were save2d
end,if nf>1, fprintf('\n'); end

