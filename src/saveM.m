function saveM(M,outdir)

if nargin~=2 || ischar(M) || ~ischar(outdir), error('usage:  saveM(M,outdir)'); end
if ~isstruct(M), error('M should be a struct'); end

if exist(outdir,'file') && ~exist(outdir,'dir')
  error('%s already exists and is not a directory',outdir);
end

if exist(outdir,'dir'), rmdir(outdir,'s'); end
ede(outdir);

f = fieldnames(M);
for i=1:length(f)
 try 
  tmp = getfield(M,f{i});

  % don't try to save FixedWidthBinary file handles
  cls = class(tmp);
  if strcmp(cls,'org.broadinstitute.cga.tools.seq.FixedWidthBinary')
    continue
  end

  if isstruct(M.(f{i})),
    save2(tmp,[outdir '/' f{i}]);
  else
    save([outdir '/' f{i} '.mat'],'tmp');
  end
 catch me
   disp(me)
   disp(me.message)
   fprintf('WARNING: failed to save %s\n',f{i});
 end
end

