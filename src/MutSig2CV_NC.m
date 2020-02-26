function [M P] = MutSig2CV_NC(M,outdir,P)
% M = MutSig2CV_v3_test3_wrapper(maf)
% M = MutSig2CV_v3_test3_wrapper(maf,P)
% [M P] = MutSig2CV_v3_test3_wrapper(maf)
% [M P] = MutSig2CV_v3_test3_wrapper(maf,P)
% G = MutSig2CV_v3_test3_wrapper(maf,outdir)
% G = MutSig2CV_v3_test3_wrapper(maf,outdir,P)
% G = MutSig2CV_v3_test3_wrapper(maf,outdir)
% G = MutSig2CV_v3_test3_wrapper(maf,outdir,P)

if nargin>3, error('too many input arguments'); end
if nargin==0, error('input required'); end
if nargout>2, error('too many output arguments'); end

if isstruct(M) && isfield(M,'mut') && isfield(M,'geneidx')
  % It's already loaded: just need to RUN
  if ~exist('outdir','var'), outdir=[]; end
  if nargin==2 && isstruct(outdir)
    P=outdir;
    outdir=[];
  end
  if ~exist('P','var'), P=[]; end
  M = MutSig2CV_v5_run(M,outdir,P);
  return
elseif isstruct(M) && isfield(M,'chr') && isfield(M,'patient')
  % maf is loaded, just need to convert to M
  maf = M;
elseif ischar(M)
  maf = M;
  demand_file(maf);
else
  error('invalid input: first argument should be maf name or loaded maf or loaded M struct');
end

load_only = false;
if nargin==1
  load_only = true;
  outdir=[];
  P=[];
elseif nargin==2 && (isstruct(outdir)||isempty(outdir))
  load_only = true;
  P = outdir;
  outdir = [];
end

if ~exist('outdir','var') || isempty(outdir) || strcmpi('none',outdir), outdir = ['/tmp/mutsig_temp' num2str(round(now*1e6+1000*rand))]; fprintf('Writing to temporary outdir %s\n',outdir); end
ede(outdir);

if exist('P','var') && ~isempty(P)
  if ischar(P)
    P = process_params_file([],P);
  elseif isstruct(P)
     % do nothing
  else
    error('unknown format for P');
  end
end

if ~exist('P','var'), P=[]; end

% FIXED parameters
P = impose_default_value(P,'context_and_effect_fwb_file',    'ref/c65e15/all.fwb');
P = impose_default_value(P,'context_and_effect_categs_file', 'ref/c65e15/categs.txt');
P = impose_default_value(P,'effect_table',                   'ref/cons_effect_table.txt');
P = impose_default_value(P,'FixedWidthBinary_jar_file',      'src/FixedWidthBinary.jar');
P = impose_default_value(P,'enforce_target_list_footprint_flank_size', 3000);

P = impose_default_value(P,'load_only',false);
if P.load_only, load_only = true; end

% REQUIRED parameters with no defaults
demand_fields(P, {'coverage_models_mat_file' 'covariates_file'})

params_file = [outdir '/params.txt'];
write_params_file(P,params_file);

args = {...
    maf ...
    outdir ...
    P.effect_table ...
    P.coverage_models_mat_file ...
    P.context_and_effect_fwb_file ...
    P.context_and_effect_categs_file ...
    P.covariates_file ...
    P.FixedWidthBinary_jar_file ...
    params_file ...
};

if load_only
  fprintf('WILL LOAD DATA ONLY AND RETURN.\n');
else
  fprintf('WILL LOAD+RUN.\n');
end

% LOAD 
[M P] = MutSig2CV_v5_load(args{:});
if load_only, return; end

% SAVE quick-load version for later debugging
%if ~grepm('\.M$',{maf})
  fname = [outdir '/input_data.M'];
  fprintf('Saving quick-load version of input data as %s       ',fname);
  tic;saveM(M,fname);toc
%end

% RUN 
M = MutSig2CV_v5_run(M,outdir,P);





