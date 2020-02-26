function G = MutSig2CV_v5_run(M,outdir,P)

if nargin>3, error('extraneous input'); end

%if ~exist('M','var') || ~isstruct(M) || isfield(M,'ng')
%  error('M should be from mload.  Maybe you meant mwrap?\n');
%end

if nargin==2 & (isstruct(outdir)||isempty(outdir))
  P=outdir;
  outdir=[];
end

% ensure output directory OK
if ~exist('outdir','var') || isempty(outdir) || strcmpi('none',outdir)
  outdir = ['/tmp/mutsig_temp' num2str(round(now*1e6+1000*rand))];
  fprintf('Writing to temporary outdir %s\n',outdir);
end
ede(outdir);

if ~exist('P','var'), P=[]; end

% MutSig RUN parameters
P = impose_default_value(P,'scatter_jobcount',1);
P = impose_default_value(P,'scatter_jobno',1);
P = impose_default_value(P,'genes_to_analyze',{});
P = impose_default_value(P,'patients_to_analyze',{});
P = impose_default_value(P,'ttypes_to_analyze',{});
P = impose_default_value(P,'genes_to_exclude',{});
P = impose_default_value(P,'patients_to_exclude',{});
P = impose_default_value(P,'ttypes_to_exclude',{});
P = impose_default_value(P,'exclude_indels',false);
P = impose_default_value(P,'exclude_noncoding',true);
if ~P.exclude_noncoding, fprintf('NOTE: including noncoding in background estimates.\n'); end
P = impose_default_value(P,'use_NMF',false);
P = impose_default_value(P,'use_new_approach',true);
P = impose_default_value(P,'use_new_prior',false);
P = impose_default_value(P,'use_simple_prior',false);
P = impose_default_value(P,'use_dense_Fg_domain',false);
P = impose_default_value(P,'use_joint_Fg_mupc_estimation',false);
P = impose_default_value(P,'n_iter_in_joint_est',2);
P = impose_default_value(P,'use_prior_in_joint_estimation',false);
P = impose_default_value(P,'use_joint_est_as_Fg_prior',false);
P = impose_default_value(P,'use_new_gg_method',false);
P = impose_default_value(P,'use_v3_test4_method',false);
P = impose_default_value(P,'num_sampling_iterations',100);
P = impose_default_value(P,'Fg_distribution_type','beta');
if ~ismember(P.Fg_distribution_type,{'beta','normal','t'}), error('invalid P.Fg_distribution_type'); end
if strcmp(P.Fg_distribution_type,'t') && ~isfield(P,'t_distribution_v'), fprintf('WARNING: using default v=1 for t-distribution!'); end
P = impose_default_value(P,'t_distribution_v',1);
P = impose_default_value(P,'use_weighted_bagels',false);
P = impose_default_value(P,'include_noncoding_in_bagel_stoptest',~P.exclude_noncoding);
P = impose_default_value(P,'include_nonsilent_in_bagel_stoptest',false);
P = impose_default_value(P,'include_noncoding_in_Fg_estimate',~P.exclude_noncoding);
P = impose_default_value(P,'include_nonsilent_in_Fg_estimate',false);
P = impose_default_value(P,'analyze_silent_as_signal',false);
P = impose_default_value(P,'max_X_in_Fg_estimate',inf);
P = impose_default_value(P,'mu_0_adjustment',0);
P = impose_default_value(P,'skip_bagels',false);
P = impose_default_value(P,'min_neighbors',0);
P = impose_default_value(P,'max_neighbors',5000);
P = impose_default_value(P,'min_mutations_in_bagel',0);
P = impose_default_value(P,'qual_min',0.1);
P = impose_default_value(P,'min_territory_ratio',0);

%P = impose_default_value(P,'indel_min_neighbors',P.min_neighbors);                       % OBSOLETE
%P = impose_default_value(P,'indel_max_neighbors',P.max_neighbors);                       % OBSOLETE
%P = impose_default_value(P,'indel_qual_min',P.qual_min);                                 % OBSOLETE
%P = impose_default_value(P,'indel_min_territory_ratio',P.min_territory_ratio);           % OBSOLETE

P = impose_default_value(P,'use_pCV_gaussian_approximation',false);
P = impose_default_value(P,'pCV_gaussian_approximation_min_npat_exp',5);
P = impose_default_value(P,'base_min_effect_size',1.01);
if P.base_min_effect_size<1, error('P.base_min_effect_size is multiplicative, must be >=1'); end
P = impose_default_value(P,'penalty_per_strike',1.10);
if P.penalty_per_strike<1, error('P.penalty_per_strike is multiplicative, must be >=1'); end
P = impose_default_value(P,'impute_full_cov_when_promotes_significance',false);

P = impose_default_value(P,'coding_BMR_threshold_for_pCL_exclusion',30e-6);
if ~isempty(P.coding_BMR_threshold_for_pCL_exclusion) && (P.coding_BMR_threshold_for_pCL_exclusion<1e-9 || P.coding_BMR_threshold_for_pCL_exclusion>1e-2), error('invalid P.coding_BMR_threshold_for_pCL_exclusion'); end 
P = impose_default_value(P,'remove_ULTRA_mutations_early',false);
P = impose_default_value(P,'restrict_to_one_mutation_per_patient',true);
P = impose_default_value(P,'radius_to_impute_coverage_around_mutations',10);
if isfield(P,'permutations_min_effect_size'), error('No longer supported: P.permutations_min_effect_size.   Please use P.pCL_min_effect_size and P.pFN_min_effect_size.'); end
P = impose_default_value(P,'pCL_min_effect_size',1.01);
if P.pCL_min_effect_size<1, error('pCL_min_effect_size must be at least 1.00'); end
P = impose_default_value(P,'pFN_min_effect_size',1.001);
if P.pFN_min_effect_size<1, error('pFN_min_effect_size must be at least 1.00'); end
P = impose_default_value(P,'pCF_min_effect_size',1.01);
if P.pCF_min_effect_size<1, error('CF_min_effect_size must be at least 1.00'); end
P = impose_default_value(P,'max_coverage_bins',10);
P = impose_default_value(P,'clustering_metric',204);
P = impose_default_value(P,'randseed',6789);
P = impose_default_value(P,'skip_pCF',true);
P = impose_default_value(P,'skip_pCV_convolutions',false);
P = impose_default_value(P,'skip_permutations',true);
P = impose_default_value(P,'maxperm',1e6); %to speedup initial screens; need to reset to 1e6 later.
if P.maxperm==0, P.skip_permutations = true; end
P = impose_default_value(P,'theta',1);
P = impose_default_value(P,'keyboard_before_begin',false);
P = impose_default_value(P,'skip_output_files',false);
P = impose_default_value(P,'skip_output_maf',false);
P = impose_default_value(P,'skip_verbose_output_mat',false);
P = impose_default_value(P,'memory_audit',false);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PREPARE FOR RUN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%memaudit

fprintf('PREPARING FOR RUN\n');

% GENES to analyze
if isempty(P.genes_to_analyze) && isempty(P.genes_to_exclude) && ~(P.scatter_jobcount>=2 && P.scatter_jobno>=1)
  genes_to_analyze = 1:M.ng;
else
  if (P.scatter_jobcount>=2 && P.scatter_jobno>=1) && (~isempty(P.genes_to_analyze)||~isempty(P.genes_to_exclude))
    error('please specify EITHER scatter jobno/count OR genes to analyze/exclude');
  elseif (P.scatter_jobcount>=2 && P.scatter_jobno>=1)
    if ~(P.scatter_jobcount>=2 && P.scatter_jobno>=1 && P.scatter_jobno<=P.scatter_jobcount), error('Please specify valid jobcount and jobno'); end
    [fg lg] = calculate_work_split(P.scatter_jobno, P.scatter_jobcount, M.ng);
    fprintf('SCATTER JOB %d/%d    will analyze genes %d-%d\n',P.scatter_jobno, P.scatter_jobcount, fg, lg);
    if ~grepm('/part\d+$',{outdir})
      outdir = [outdir '/part' num2str(P.scatter_jobno)];
    else
      tmp = parse({outdir},'/part(\d+)$',{'no'},1);
      if tmp.no~=P.scatter_jobno, error('outdir should end with /part%d',P.scatter_jobno); end
    end
    genes_to_analyze = fg:lg;
  elseif ~isempty(P.genes_to_analyze) && ~isempty(P.genes_to_exclude)
    error('please specify EITHER P.genes_to_analyze OR P.genes_to_exclude');
  elseif ~isempty(P.genes_to_analyze)
    if ischar(P.genes_to_analyze), P.genes_to_analyze = {P.genes_to_analyze}; end
    if isnumeric(P.genes_to_analyze)
      genes_to_analyze = P.genes_to_analyze;
    elseif iscellstr(P.genes_to_analyze)
      genes_to_analyze = listmap(P.genes_to_analyze,M.gene.name);
    else
      error('invalid P.genes_to_analyze');
    end
  elseif ~isempty(P.genes_to_exclude)
    if ischar(P.genes_to_exclude), P.genes_to_exclude = {P.genes_to_exclude}; end
    if isnumeric(P.genes_to_exclude)
      genes_to_analyze = setdiff(1:M.ng,P.genes_to_exclude);
    elseif iscellstr(P.genes_to_exclude)
      genes_to_analyze = find(~ismember(M.gene.name,P.genes_to_exclude));
    else
      error('invalid P.genes_to_exclude');
    end
  end
end
genes_to_analyze = as_row(genes_to_analyze(~isnan(genes_to_analyze)&genes_to_analyze>=1&genes_to_analyze<=M.ng));
if isempty(genes_to_analyze) error('invalid P.genes_to_analyze/exclude'); end

% OUTPUT DIRECTORY (note: may have been edited to ".../partNNN" above)
ede(outdir);
if exist([outdir '/results.mat'],'file')
  fprintf('WARNING: results already exist in this directory!\n');
  fprintf('Pausing 10 seconds before continuing.\n');
  pause(10);
end
save([outdir '/full_params.mat'],'P');
logfile = fopen([outdir '/log.txt'],'wt');
msgfile = fopen([outdir '/msg.txt'],'wt');

% DATA SUBSETTING

%if P.exclude_indels || P.exclude_noncoding || ~isempty(P.patients_to_analyze) || ~isempty(P.ttypes_to_analyze) ||...
%  ~isempty(P.patients_to_exclude) || ~isempty(P.ttypes_to_exclude)
%  if (~isempty(P.patients_to_analyze))+(~isempty(P.ttypes_to_analyze))+(~isempty(P.patients_to_exclude))+(~isempty(P.ttypes_to_exclude))>1
%    error('Please specify ONLY ONE of P.patients_to_analyze, P.ttypes_to_analyze, P.patients_to_exclude, P.ttypes_to_exclude');
%  end
%  if ~isempty(P.patients_to_analyze) || ~isempty(P.ttypes_to_analyze) || ~isempty(P.patients_to_exclude) || ~isempty(P.ttypes_to_exclude)
%    if ~isempty(P.ttypes_to_analyze)
%      if ~isfield(M.pat,'ttype'), error('Missing ttype information'); end
%      if ischar(P.ttypes_to_analyze), P.ttypes_to_analyze = {P.ttypes_to_analyze}; end
%      if iscellstr(P.ttypes_to_analyze)
%        if ismember(P.ttypes_to_analyze,'PANCAN')
%          patients_to_analyze = 1:M.np;
%        else
%          patients_to_analyze = find(ismember(M.pat.ttype,P.ttypes_to_analyze));
%        end
%      else
%        error('invalid P.ttypes_to_analyze');
%      end
%    elseif ~isempty(P.ttypes_to_exclude)
%      if ~isfield(M.pat,'ttype'), error('Missing ttype information'); end
%      if ischar(P.ttypes_to_exclude), P.ttypes_to_exclude = {P.ttypes_to_exclude}; end
%      if iscellstr(P.ttypes_to_exclude)
%        if ismember(P.ttypes_to_exclude,'PANCAN')
%          error('can''t exclude PANCAN');
%        else
%          patients_to_analyze = find(~ismember(M.pat.ttype,P.ttypes_to_exclude));
%        end
%      else
%        error('invalid P.ttypes_to_exclude');
%      end
%    elseif ~isempty(P.patients_to_analyze)
%      if ischar(P.patients_to_analyze), P.patients_to_analyze = {P.patients_to_analyze}; end
%      if isnumeric(P.patients_to_analyze)
%        patients_to_analyze = P.patients_to_analyze;
%      elseif iscellstr(P.patients_to_analyze)
%        patients_to_analyze = listmap(P.patients_to_analyze,M.pat.name);
%      else
%        error('invalid P.patients_to_analyze');
%      end
%    elseif ~isempty(P.patients_to_exclude)
%      if ischar(P.patients_to_exclude), P.patients_to_exclude = {P.patients_to_exclude}; end
%      if isnumeric(P.patients_to_exclude)
%        patients_to_analyze = setdiff(1:M.np,P.patients_to_exclude);
%      elseif iscellstr(P.patients_to_exclude)
%        patients_to_analyze = find(~ismember(M.pat.name,P.patients_to_exclude));
%      else
%        error('invalid P.patients_to_analyze');
%      end
%    end
%    fprintf('Restricting to analysis of %d patients.\n', length(patients_to_analyze));
%    M.pat.idx_orig = as_column(1:M.np);
%    M.pat = reorder_struct(M.pat,patients_to_analyze);
%    M.np = slength(M.pat);
%    M.mut = reorder_struct(M.mut,ismember(M.mut.pat_idx,patients_to_analyze));
%    M.mut.pat_idx = mapacross(M.mut.pat_idx,M.pat.idx_orig,1:M.np);
%    M.pat = rmfield(M.pat,'idx_orig');
%  end
%  if P.exclude_indels
%    fprintf('Excluding indels.\n');
%    M.mut = reorder_struct_exclude(M.mut,ismember(M.mut.effect_idx,find(grepmi('indel',M.effect.name))));
%  end
%  fprintf('Recomputing geneidx lookup table...\n');
%  M.mut = sort_struct(M.mut,'gene_idx');
%  [u ui uj] = unique(M.mut.gene_idx,'first');
%  h = histc(uj,1:length(u));
%  M.geneidx = cell(M.ng,1); for i=1:length(u), M.geneidx{u(i)} = as_column(ui(i):ui(i)+h(i)-1);end
%end

% (re-)open FWB tracks
if ~P.skip_permutations,
  M.FWB.context_and_effect = org.broadinstitute.cga.tools.seq.FixedWidthBinary(M.FWB.context_and_effect_file);
end

fprintf('Processing covariates...\n');
% convert covariate raw values to Z-scores
M.Z = nan(M.ng,M.nv);
for vi=1:M.nv
  missing = isnan(M.V(:,vi)) | isinf(M.V(:,vi));
  mn = mean(M.V(~missing,vi));
  sd = std(M.V(~missing,vi),0);  % second parameter=0 means confirm default behavior of normalize by (N-1) not (N)
  M.Z(~missing,vi) = (M.V(~missing,vi)-mn)./sd;
end

% compute number of "strikes"
if all(isfield(M.gene,{'log_exprmax','rt','hiC','paz'}))
  % legacy method
  M.gene.nstrikes = (M.gene.log_exprmax<4.5)+(M.gene.log_exprmax<4.0)+(M.gene.rt>600)+(M.gene.rt>800)+...
      (M.gene.hiC<-0.02)+(M.gene.hiC<-0.01)+(M.gene.paz>0.2)+(M.gene.paz>0.3);
else
  fprintf('*********************************************************************************************\n');
  fprintf('WARNING: these are not the standard covariates.  Using UNTESTED METHOD of computing nstrikes.\n');
  fprintf('*********************************************************************************************\n');
  % approximate general method for future release, NOT TESTED YET
  % (needs to be replaced with a bootstrapping method)
  M.gene.nstrikes = sum(M.Z>1.5,2) + sum(M.Z>2,2);
end

% compute per-gene minimum effect size
M.gene.min_effect_size = P.base_min_effect_size*P.penalty_per_strike.^(M.gene.nstrikes);

% PERMUTATIONS preparation
if ~P.skip_permutations
  rand('twister',P.randseed);   % initialize random number generator
  M.mut.trackpos = nan(slength(M.mut),1);    % allocate field for "position along track" that will be filled out during permutations
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BACKGROUND MODEL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TOTAL COUNTS

%memaudit

fprintf('Computing counts and rates...\n');

% increase territory to the max coverage represented (this prevents a class of edge-case problems later)
M.cov.CE_terr = ceil(max([M.cov.CE_terr M.cov.CE_cov],[],2));

%note this is never used except for output purposes.
M.gene.codelen = round(sum(sum(M.cov.CE_terr,3),2)/3);

%M.gene.[Nn].* are only used for bagel calculation purposes.
M.gene.Nssnv = round(squeeze(sum(M.cov.CE_cov, 2))*M.np);
M.gene.Nind = sum(M.gene.Nssnv, 2);

%M.gene.Nncd=0;
%M.gene.Nsyn=0;
%M.gene.Nmis=0;
%M.gene.Nnon=0;
%M.gene.Nspl=0;
%for cov_idx=1:M.cov.ntype
%  pidx = find(M.pat.cov_idx==cov_idx);
%  M.gene.Nncd = M.gene.Nncd + sum(M.cov.gene_effect_cov(:,cov_idx,:,1),3)*length(pidx);
%  M.gene.Nsyn = M.gene.Nsyn + sum(M.cov.gene_effect_cov(:,cov_idx,:,2),3)*length(pidx);
%  M.gene.Nmis = M.gene.Nmis + sum(M.cov.gene_effect_cov(:,cov_idx,:,3),3)*length(pidx);
%  M.gene.Nnon = M.gene.Nnon + sum(M.cov.gene_effect_cov(:,cov_idx,:,4),3)*length(pidx);
%  M.gene.Nspl = M.gene.Nspl + sum(M.cov.gene_effect_cov(:,cov_idx,:,5),3)*length(pidx);
%end
%M.gene.Nind = M.gene.Nncd + M.gene.Nsyn + M.gene.Nmis + M.gene.Nnon + M.gene.Nspl;
%
%M.gene.Nncd = round(M.gene.Nncd);
%M.gene.Nsyn = round(M.gene.Nsyn);
%M.gene.Nmis = round(M.gene.Nmis);
%M.gene.Nnon = round(M.gene.Nnon);
%M.gene.Nspl = round(M.gene.Nspl);
%M.gene.Nind = round(M.gene.Nind);

n_all = full(hist2d_sparse(M.mut.gene_idx, M.mut.effect_idx, 1, M.ng, 1, slength(M.effect)));
M.gene.nssnv = n_all(:, ~M.effect.is_indel);
M.gene.nind = n_all(:, M.effect.is_indel);

%ncd = (M.mut.effect_idx==find(strcmp('ncd',M.effect.name)));
%syn = (M.mut.effect_idx==find(strcmp('syn',M.effect.name)));
%mis = (M.mut.effect_idx==find(strcmp('mis',M.effect.name)));
%non = (M.mut.effect_idx==find(strcmp('non',M.effect.name)));
%spl = (M.mut.effect_idx==find(strcmp('spl',M.effect.name)));
%ind = (ismember(M.mut.effect_idx,find(grepmi('indel_cod|indel_spl',M.effect.name))));
%
%M.gene.nncd = as_column(histc(M.mut.gene_idx(ncd),1:M.ng));
%M.gene.nsyn = as_column(histc(M.mut.gene_idx(syn),1:M.ng));
%M.gene.nmis = as_column(histc(M.mut.gene_idx(mis),1:M.ng));
%M.gene.nnon = as_column(histc(M.mut.gene_idx(non),1:M.ng));
%M.gene.nspl = as_column(histc(M.mut.gene_idx(spl),1:M.ng));
%M.gene.nind = as_column(histc(M.mut.gene_idx(ind),1:M.ng));

% distribution of gene rates
% note that these aren't used anywhere with current Fg estimation scheme.
%M.gene.rate_sil = M.gene.nsyn ./ M.gene.Nsyn;
%M.gene.rate_tot = (M.gene.nsyn+M.gene.nmis+M.gene.nnon+M.gene.nspl) ./ (M.gene.Nsyn+M.gene.Nmis+M.gene.Nnon+M.gene.Nspl);
%mu_0 = median(M.gene.rate_tot)*(1+P.mu_0_adjustment);
%if mu_0<0.5e-6, mu_0 = mean(M.gene.rate_tot); end

%memaudit

% PATIENT RATES 

M.pat.N_tot = ones(M.np, 1)*fullsum(M.cov.CE_cov);
M.pat.N_tot = round(M.pat.N_tot/3);
M.pat.nbg_tot = as_column(histc(M.mut.pat_idx(ismember(M.mut.effect_idx, find(M.effect.is_background))),1:M.np));
M.pat.nsig_tot = as_column(histc(M.mut.pat_idx(ismember(M.mut.effect_idx, find(~M.effect.is_background))),1:M.np));
M.pat.n_tot = as_column(histc(M.mut.pat_idx,1:M.np));
M.pat.rate_bg = M.pat.nbg_tot./M.pat.N_tot;  % (not used anywhere)
M.pat.rate_sig = M.pat.nsig_tot./M.pat.N_tot;  % (not used anywhere)
M.pat.rate_tot = M.pat.n_tot./M.pat.N_tot; % used to exclude high mutation rate patients from pCL/FN
M.pat.log_rate_tot = max(-9,log10(M.pat.rate_tot)); % was used in joint mu_pc/Fg estimate; not used anymore.

% also compute table of (patient,category)-specific rates
M.pat.n_c = hist2d_fast(M.mut.pat_idx,M.mut.categ_idx,1,M.np,1,M.ncat);
midx = ismember(M.mut.effect_idx,find(M.effect.is_indel));
M.pat.n_ind = as_column(histc(M.mut.pat_idx(midx),1:M.np));
%M.pat.N_c_2 = repmat(squeeze(sum(sum(M.cov.gene_effect_terr(:,1,:,2:5),1),4))',M.np,1);
M.pat.N_c_2 = repmat(squeeze(sum(sum(M.cov.CE_cov,1),3)),M.np,1);
M.pat.N_ind = sum(M.pat.N_c_2,2);
if P.use_NMF, M.pat.n_c = M.pat.n_c + 0.0001; end    % To run nmf with patients with no coding mutations
M.pat.mu_c = M.pat.n_c ./ M.pat.N_c_2;
M.pat.mu_ind = M.pat.n_ind ./ M.pat.N_ind;

% NMF Smoothing of channels
if P.use_NMF
  fprintf('NMF Smoothing using k=%d\n',P.use_NMF);
  [h1 w1] = nmf(M.pat.mu_c,P.use_NMF,0);
  M.pat.mu_c = h1*w1;
end

% EXCLUSION OF HIGH-MUTATION-RATE PATIENTS FROM pCL
M.pat.exclude_from_pCL = false(slength(M.pat),1);
if ~isempty(P.coding_BMR_threshold_for_pCL_exclusion)
  M.pat.exclude_from_pCL = (M.pat.rate_tot>P.coding_BMR_threshold_for_pCL_exclusion);
  fprintf('Identified %d/%d hypermutated patients (>%0.2d/Mb mutations) to be excluded from pCL\n',sum(M.pat.exclude_from_pCL),slength(M.pat),1e6*P.coding_BMR_threshold_for_pCL_exclusion);
end

%memaudit

% PRIORS 

if P.use_dense_Fg_domain
  Fg_domain = as_column(10.^[-1:0.05:1]);
else
  Fg_domain = as_column(10.^[-1 -0.5 -0.2 0 0.2 0.5 1]);
end
nfg = length(Fg_domain);

if (P.use_simple_prior+P.use_new_prior+P.use_joint_est_as_Fg_prior)>1
  error('please use only one');
end

if P.use_simple_prior || P.use_new_prior || P.use_joint_est_as_Fg_prior
  if ~P.use_new_approach, error('Priors only implemented within "new approach"'); end
end
  
if P.use_simple_prior
  fprintf('use_simple_prior\n');
  if P.use_dense_Fg_domain
    p = normpdf(log10(Fg_domain),-0.1,0.21);
    Fg_prior = as_row(p/sum(p));
  else
    Fg_prior = [0.0095    0.0702    0.3286    0.3402    0.2120    0.0394    0.0003];
  end
end

if P.use_new_prior
  fprintf('use_new_prior: Computing P(Fg|data) for each gene...\n');
  all_Fg_prob = nan(M.ng,nfg);
  for g=genes_to_analyze, if ~mod(g,1000), fprintf('%d ',g); end
    gene_midx = M.geneidx{g};
    midx0 = gene_midx(M.mut.effect_idx(gene_midx)==find(strcmp('syn',M.effect.name)));
    gene_nsyn = hist2d_fast(M.mut.pat_idx(midx0),M.mut.categ_idx(midx0),1,M.np,1,M.ncat);
    gene_Nsyn = zeros(M.np,M.ncat);
    for c=1:M.ncat
      gene_Nsyn(:,c) = squeeze(M.cov.gene_effect_cov(g,M.pat.cov_idx,c,2));
    end
    gene_nsyn = round(gene_nsyn);
    gene_Nsyn = round(gene_Nsyn);
    gene_Nsyn = max(gene_Nsyn,10*gene_nsyn);     % (to kill N==0 cases)
    lp = nan(nfg,1);
    for fgi=1:nfg
      mu = Fg_domain(fgi) * M.pat.mu_c;
      lp(fgi) = nansum(log(binopdf(gene_nsyn(:),gene_Nsyn(:),mu(:))));
    end
    lp = lp - max(lp);
    p = exp(lp);
    p = p/sum(p);
    all_Fg_prob(g,:) = p;
  end, fprintf('\n');
  % BUILD A PRIOR FROM THESE DISTRIBUTIONS
  y = all_Fg_prob;
  covar_to_use = 'exprmax';   % rt  -exprmax  -hiC
  gene_window_size = 2000;
  [cvar_sorted ord] = sort(M.gene.(covar_to_use));
  y=y(ord,:);
  z = nan(size(y));
  st=1; en=gene_window_size; mid=round((st+en)/2);
  for i=1:M.ng
    z(i,:) = nansum(y(st:en,:),1);
    if i>mid & en<M.ng, st=st+1; mid=mid+1; en=en+1; end
  end
  z = bsxfun(@rdivide,z,sum(z,2));
  M.Fg_prior = nan(M.ng,nfg);
  M.Fg_prior(ord,:) = z;
end

%memaudit

% BAGELS

if P.max_neighbors==0
  P.skip_bagels = true;
end

M.gene.ntest = M.gene.nssnv(:, M.effect.is_background);
M.gene.Ntest = M.gene.Nssnv(:, M.effect.is_background);

% will discover two sets of bagels:
%      main_bagel
%      indel_bagel (OBSOLETE)
%n_bageltypes = 2;
%min_neighbors = [P.min_neighbors P.indel_min_neighbors];
%max_neighbors = [P.max_neighbors P.indel_max_neighbors];
%min_territory_ratio = [P.min_territory_ratio P.indel_min_territory_ratio];
%qual_min = [P.qual_min P.indel_qual_min];

n_bageltypes = 1;
min_neighbors = [P.min_neighbors];
max_neighbors = [P.max_neighbors];
min_territory_ratio = [P.min_territory_ratio];
qual_min = [P.qual_min];

max_min_neighbors = max(min_neighbors);
max_max_neighbors = max(max_neighbors);
max_min_territory_ratio = max(min_territory_ratio);
min_qual_min = min(qual_min);

M.gene.nbagel = zeros(M.ng,1);
M.gene.Nbagel = zeros(M.ng,1);
M.gene.nnei = zeros(M.ng,1);
M.gene.nnei2 = M.gene.nnei;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% RUN 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%memaudit

z = nan(M.ng,1);
z2 = nan(M.ng,nfg);
M.gene.npat_exp = z;
M.gene.npat = z; M.gene.nsite = z;
M.gene.x = z; M.gene.X = z;
M.gene.Fg_post = z2;
M.gene.S = nan(M.ng,M.np);
M.gene.pCVmin = z; M.gene.pCVmid = z; M.gene.pCVmax = z;
M.gene.nmut = z;
M.gene.sCL = z; M.gene.sFN = z;

M.gene.pFN_Gaussian_m = z;
M.gene.pFN_Gaussian_s = z;
M.gene.pFN_ocons = z;

M.gene.pCV_Gaussian = z2;
M.gene.pCV_Gaussian_m = z2;
M.gene.pCV_Gaussian_s = z2;
M.gene.pCV_score_obs = z2;
M.gene.pCV_mefsz = z2;

M.gene.nperm = z;
M.gene.pCL = z; M.gene.pFN = z; M.gene.pCLFN = z;
M.gene.pCL2 = z; M.gene.pFN2 = z; M.gene.pCLFN2 = z;
M.gene.pCF = z;

%%%%%%%%%%%%%%%%%%%
% pCF
% --> use nonsilent coding mutations only
% --> use raw AF (i_tumor_f)
%%%%%%%%%%%%%%%%%%%

if P.skip_pCF
  fprintf('Skipping calculation of pCF.\n');
else

  if ~isfield(M.mut,'i_tumor_f')
    fprintf('Can''t compute pCF because i_tumor_f field not present.\n');
  else
    m1 = reorder_struct(M.mut,mis|non|spl|ind);
    m1 = make_numeric(m1,'i_tumor_f');
    m1 = reorder_struct(m1,m1.i_tumor_f>=0 & m1.i_tumor_f<=1);
    if slength(m1)==0
      fprintf('Can''t compute pCF because i_tumor_f field is empty/invalid.\n');
    else
      try
        fprintf('Computing pCF...\n');
        
% OLD STYLE
%        m1.pctaf = 100 * m1.i_tumor_f;
%        PP=[]; PP.muts_to_use = 'both'; PP.coding_only = false; % (don't try to do any more filtering)
%        PP.min_effect_size = P.pCF_min_effect_size;
%        tmp=[]; [tmp.p tmp.q tmp.gene tmp.chr] = calc_af_score(m1,PP);
%        M.gene.pCF = mapacross(M.gene.name,tmp.gene,tmp.p);

% NEW STYLE
        PP=[]; PP.min_effect_size = P.pCF_min_effect_size;
        M.gene.pCF = calc_pCF(m1,M.ng,M.np,PP);

        fprintf('pCF successfully computed!\n');
      catch me
        fprintf('ERROR WHILE CALCULATING pCF\n');
        disp(me.message)
        disp(me)
        for i=1:length(me.stack), disp(me.stack(i)); end
      end
    end      
  end

end




%%%%%%%%%%%%%%%%%%%
% pCV pCL pFN
%%%%%%%%%%%%%%%%%%%

fprintf('\nProcessing each gene...\n');

fprintf('now         eta                   # gene            nmut   nperm  CV                    CL                    FN                     [ qq plot                       ]  time_usage\n');

if P.keyboard_before_begin, keyboard; end

qqhist = zeros(1,6);  % for 'live' QQ statistics

% for time profiling of the gene loop
tt=tic; timeused=zeros(1,6);

M.gene.shortname = regexprep(M.gene.name, '.*gencode::(.*)::.*', '$1');

%genes_to_analyze = find(ismember(M.gene.shortname, {'TP53' 'PIK3CA' 'CBFB' 'MAP3K1' 'PTEN' 'GATA3' 'AKT1' 'SF3B1' 'CDH1' 'CABP2' 'XBP1' 'PTRHD1'}))';
%genes_to_analyze = find(ismember(M.gene.shortname, {'RAB40C'}))';

%genes_to_analyze = find(strcmp(M.gene.name, 'enhancers::chr6:142705600-142706400::NA::NA'));

%genes_to_analyze = [417 973];

%genes_to_analyze = 7094;

num_genes_analyzed=0;
for g=genes_to_analyze
 %if mod(g,1000)==1, memaudit; end

 gene_status_reported=false;
 pCVmax = 1; pclust=1; pcons=1; nperm=0; nm=0;
 try

 if P.use_v3_test4_method

   pCVmax = avg_p(g); pCVmid = avg_p_mid(g); pCVmin = avg_p_next(g);

 else

 if ~P.skip_bagels

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % BAGEL FINDING
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ttt=tic;
  % calculate distances from this gene
  df2 = bsxfun(@minus,M.Z,M.Z(g,:)).^2; dist2 = nansum(df2,2)./sum(~isnan(df2),2); [tmp,ord] = sort(dist2); ord = [g;ord(ord~=g)];
  per_type_nnei = zeros(n_bageltypes,1); per_type_bagel = repmat({nan(1,max_max_neighbors)},n_bageltypes,1);
  per_type_bagel_done = false(n_bageltypes,1);
  % expand bagel outward until quality falls below qual_min
  nfit=0; Nfit=0; nbagel=0; Nbagel=0;
  for ni=0:max_max_neighbors, gidx = ord(ni+1);
    ngene = M.gene.ntest(gidx); Ngene = M.gene.Ntest(gidx);
    if ni==0, ngene0=ngene; Ngene0=Ngene; end
    nfit=nfit+ngene; Nfit=Nfit+Ngene;
    % compare the gene being added to the central gene
    qual = 2*hyge2cdf(ngene,Ngene,ngene0,Ngene0); if qual>1, qual = 2-qual; end
    territory_ratio = Nfit/Ngene0;
    % stopping criterion: stop if this gene would drop quality below qual_min
    if ((ni-1)>=max_min_neighbors && territory_ratio>=max_min_territory_ratio) && qual<min_qual_min && nbagel>=P.min_mutations_in_bagel, break; end
    if ni>0   % update gene's bagel
      %M.gene.bagel(g,ni)=gidx;   % NOTE: not actually used anywhere
      for c=1:n_bageltypes    % which categories include this bagel?
        if ~per_type_bagel_done(c) && (ni<=min_neighbors(c) || territory_ratio<min_territory_ratio(c) || nbagel<P.min_mutations_in_bagel || ...
                                       (ni<=max_neighbors(c) && qual>=qual_min(c)))
          per_type_nnei(c) = per_type_nnei(c) + 1; per_type_bagel{c}(ni) = gidx;
          if c==1, nbagel=nbagel+ngene; Nbagel=Nbagel+Ngene; end
        else
          per_type_bagel_done(c) = true;
        end
      end
    end
  end % next neighborhood size
  M.gene.nbagel(g) = nbagel; M.gene.Nbagel(g) = Nbagel;
  M.gene.nnei(g) = per_type_nnei(1);   % main bagel
%  M.gene.nnei2(g) = per_type_nnei(2);  % indel bagel   (OBSOLETE)
  timeused(1)=timeused(1)+toc(ttt);

 end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % COUNTS LOOKUP
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ttt=tic;
  % signal
%  gene_midx = M.geneidx{g};
%  if P.analyze_silent_as_signal
%    midx0 = gene_midx(M.mut.effect_idx(gene_midx)==find(strcmp('syn',M.effect.name)));
%    gene_nsyn = hist2d_fast(M.mut.pat_idx(midx0),M.mut.categ_idx(midx0),1,M.np,1,M.ncat);
%  end
%  midx1 = gene_midx(M.mut.effect_idx(gene_midx)==find(strcmp('mis',M.effect.name)));
%  gene_nmis = hist2d_fast(M.mut.pat_idx(midx1),M.mut.categ_idx(midx1),1,M.np,1,M.ncat);
%  midx2 = gene_midx(M.mut.effect_idx(gene_midx)==find(strcmp('non',M.effect.name)));
%  gene_nnon = hist2d_fast(M.mut.pat_idx(midx2),M.mut.categ_idx(midx2),1,M.np,1,M.ncat);
%  midx3 = gene_midx(M.mut.effect_idx(gene_midx)==find(strcmp('spl',M.effect.name)));
%  gene_nspl = hist2d_fast(M.mut.pat_idx(midx3),M.mut.categ_idx(midx3),1,M.np,1,M.ncat);
%  midx4 = gene_midx(ismember(M.mut.effect_idx(gene_midx),find(grepm('indel_(cod|spl)',M.effect.name))));
%  if ~isempty(midx4), gene_nind = as_column(histc(M.mut.pat_idx(midx4),1:M.np)); else gene_nind = zeros(M.np,1); end

%  % background (=gene+bagel): noncoding and synonymous
%  main_sphere_nncd = zeros(M.np,M.ncat);
%  main_sphere_Nncd = zeros(M.np,M.ncat);
%  main_sphere_nsyn = zeros(M.np,M.ncat);
%  main_sphere_Nsyn = zeros(M.np,M.ncat);
%  if ~P.skip_bagels
%    bageltype = 1;
%    main_sphere = [g per_type_bagel{bageltype}(1:per_type_nnei(bageltype))];
%  else
%    main_sphere = g;
%  end
%  main_sphere_midx = cat(1,M.geneidx{main_sphere});
%  for c=1:M.ncat
%    midx = main_sphere_midx(M.mut.effect_idx(main_sphere_midx)==find(strcmp('ncd',M.effect.name)) & M.mut.categ_idx(main_sphere_midx)==c);
%    if ~isempty(midx), main_sphere_nncd(:,c) = as_column(histc(M.mut.pat_idx(midx),1:M.np)); end
%    midx = main_sphere_midx(M.mut.effect_idx(main_sphere_midx)==find(strcmp('syn',M.effect.name)) & M.mut.categ_idx(main_sphere_midx)==c);
%    if ~isempty(midx), main_sphere_nsyn(:,c) = as_column(histc(M.mut.pat_idx(midx),1:M.np)); end
%    for cov_idx=1:M.cov.ntype
%      pidx = (M.pat.cov_idx==cov_idx);
%      main_sphere_Nncd(pidx,c) = sum(M.cov.gene_effect_cov(main_sphere,cov_idx,c,1),1);
%      main_sphere_Nsyn(pidx,c) = sum(M.cov.gene_effect_cov(main_sphere,cov_idx,c,2),1);
%    end
%  end

  %aggregate (patient,categ)-wise background counts and territory across bagel

  main_sphere_x = zeros(M.np, M.ncat);
  main_sphere_X = zeros(M.np, M.ncat);

  if ~P.skip_bagels
    bageltype = 1;
    main_sphere = [g per_type_bagel{bageltype}(1:per_type_nnei(bageltype))];
  else
    main_sphere = g;
  end
  main_sphere_midx = cat(1,M.geneidx{main_sphere});
  for c = 1:M.ncat,
    midx = main_sphere_midx(ismember(M.mut.effect_idx(main_sphere_midx), find(M.effect.is_background)) & M.mut.categ_idx(main_sphere_midx) == c);
    if ~isempty(midx),
      main_sphere_x(:, c) = as_column(histc(M.mut.pat_idx(midx), 1:M.np));
    end

    main_sphere_X(:, c) = sum(sum(M.cov.CE_cov(main_sphere, c, find(M.effect.is_background)), 3), 1)*ones(M.np, 1);
  end

  % find # unique patients and sites (only for display in table; not used directly in calculation)
  midx = M.geneidx{g};
  midx = midx(ismember(M.mut.effect_idx(midx), find(~M.effect.is_background)));
  if isempty(midx)
    M.gene.nsite(g) = 0;
    M.gene.npat(g) = 0;
  else
    M.gene.nsite(g) = length(unique_combos(M.mut.chr(midx),M.mut.pos(midx)));
    M.gene.npat(g) = length(unique(M.mut.pat_idx(midx)));
  end

  %compute (patient,categ)-wise coverage for this gene
  %TODO: since we only have one callscheme, then there's no need for these to be repmat'd
  gene_Nbg = repmat(sum(M.cov.CE_cov(g, :, M.effect.is_background & ~M.effect.is_indel), 3), M.np, 1);
  gene_Nsig = repmat(M.cov.CE_cov(g, :, ~M.effect.is_background & ~M.effect.is_indel), M.np, 1);
  gene_Nind = repmat(sum(sum(M.cov.CE_cov(g, :, :), 3), 2), M.np, 1);

%  % gene's coverage (or territory)
%  if P.impute_full_cov_when_promotes_significance
%    % use territory
%    tmp = squeeze(M.cov.gene_effect_terr(g,1,:,2:5));
%    gene_Nsyn = repmat(tmp(:,1)',M.np,1);
%    gene_Nmis = repmat(tmp(:,2)',M.np,1);
%    gene_Nnon = repmat(tmp(:,3)',M.np,1);
%    gene_Nspl = repmat(tmp(:,4)',M.np,1);
%  else  % use coverage
%    gene_Nsyn = zeros(M.np,M.ncat);
%    gene_Nmis = zeros(M.np,M.ncat);
%    gene_Nnon = zeros(M.np,M.ncat);
%    gene_Nspl = zeros(M.np,M.ncat);
%    for c=1:M.ncat
%      gene_Nsyn(:,c) = squeeze(M.cov.gene_effect_cov(g,M.pat.cov_idx,c,2));
%      gene_Nmis(:,c) = squeeze(M.cov.gene_effect_cov(g,M.pat.cov_idx,c,3));
%      gene_Nnon(:,c) = squeeze(M.cov.gene_effect_cov(g,M.pat.cov_idx,c,4));
%      gene_Nspl(:,c) = squeeze(M.cov.gene_effect_cov(g,M.pat.cov_idx,c,5));
%    end
%  end
%  gene_Nind = sum(gene_Nsyn + gene_Nmis + gene_Nnon + gene_Nspl,2);

%  if ~P.analyze_silent_as_signal
%    nm = sum((any((gene_nmis+gene_nnon+gene_nspl)>0,2)+gene_nind)>0);  % (only used for display in status report line)
%    ndeg = 4;   % mis non spl ind
%    n_signal_ssnv = cat(3,gene_nmis,gene_nnon,gene_nspl);
%    N_signal_ssnv = cat(3,gene_Nmis,gene_Nnon,gene_Nspl);
%    n_signal_indel = gene_nind;
%    N_signal_indel = gene_Nind;
%  else % analyze silent as signal
%    nm = sum(any(gene_nsyn>0,2));
%    ndeg = 1;
%    n_signal_ssnv = gene_nsyn;
%    N_signal_ssnv = gene_Nsyn;
%  end

  gene_midx = M.geneidx{g};
  
  midx = gene_midx(ismember(M.mut.effect_idx(gene_midx), find(~M.effect.is_background & ~M.effect.is_indel)));

  nef = nnz(~M.effect.is_background & ~M.effect.is_indel);
  if isempty(midx),
    n_signal_ssnv = zeros(M.np, M.ncat, nef);
  else
    n_signal_ssnv = accumarray([M.mut.pat_idx(midx), M.mut.categ_idx(midx), M.effect.deg_ord(M.mut.effect_idx(midx))],...
			       1, [M.np M.ncat nef]);
  end 

  midx = gene_midx(ismember(M.mut.effect_idx(gene_midx), find(~M.effect.is_background & M.effect.is_indel)));

  if isempty(midx),
    n_signal_ind = zeros(M.np, 1);
  else 
    n_signal_ind = as_column(histc(M.mut.pat_idx(midx), 1:M.np));
  end

  N_signal_ssnv = max(gene_Nsig,10*n_signal_ssnv);     % (to kill N==0 cases)
  has_mutation_ssnv = squeeze(sum(n_signal_ssnv,2)>0); %marginalize over categories
  N_signal_indel = max(gene_Nind,10*n_signal_ind);  % (to kill N==0 cases)

  has_mutation_indel = n_signal_ind>0;
  has_mutation = cat(2,has_mutation_ssnv,has_mutation_indel);

  nm = sum(any(sum(n_signal_ssnv, 3), 2) + n_signal_ind);
  ndeg = size(n_signal_ssnv, 3) + 1; %NOTE: hardcoded for only one indel category.

  if P.skip_pCV_convolutions && ~P.skip_permutations && length(M.geneidx{g}) < 2,
    continue;
  end

  %note that n_signals are never used again

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Fg ESTIMATION
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if P.use_new_approach                                  % NEW APPROACH 
    x = main_sphere_x; X = main_sphere_X;
    x = round(x); X = round(X);
    X = max(X,10*x);  % (to kill N==0 cases)
    M.gene.x(g) = fullsum(x); M.gene.X(g) = fullsum(X); %only used for output purposes
    lp = nan(1,nfg);
    for fgi=1:nfg
      mu = Fg_domain(fgi) * M.pat.mu_c;
      lbp = log(binopdf(x(:),X(:),mu(:)));
      lp(fgi) = sum(lbp(~isnan(lbp)&~isinf(lbp)));
    end
    lp = lp - max(lp);
    p = exp(lp);
    if P.use_new_prior || P.use_joint_est_as_Fg_prior
      p = p .* M.Fg_prior(g,:);
    end
    if P.use_simple_prior
      p = p .* Fg_prior;
    end
    Fg_prob = p/sum(p);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif ~isfield(M.gene,'true_Fg')                       % OLD approach
    x = M.gene.nsyn(g); X = M.gene.Nsyn(g);
    if P.include_noncoding_in_Fg_estimate
      x=x+M.gene.nncd(g); X=X+M.gene.Nncd(g);
    end
    if P.include_nonsilent_in_Fg_estimate
      x=x+M.gene.nmis(g)+M.gene.nnon(g)+M.gene.nspl(g);
      X=X+M.gene.Nmis(g)+M.gene.Nnon(g)+M.gene.Nspl(g);
    end
    if P.max_neighbors>0 && M.gene.nnei(g)>0
      x=x+nbagel; X=X+Nbagel;
    end
    if X > P.max_X_in_Fg_estimate
      mu_mle = x/X;
      X = P.max_X_in_Fg_estimate;
      x = mu_mle*X;
    end
    x = round(x); X = round(X);
    M.gene.x(g) = x; M.gene.X(g) = X;
    if strcmpi(P.Fg_distribution_type,'beta')
      p = binopdf(x,X,mu_0*Fg_domain);
    elseif strcmpi(P.Fg_distribution_type,'normal')
      m = x/X;
      s = gaussian_sigma_from_x_X(x,X);
      qq = ((mu_0*Fg_domain)-m)/s;
      p = normpdf(qq,0,1);
    elseif strcmpi(P.Fg_distribution_type,'t')
      m = x/X;
      s = gaussian_sigma_from_x_X(x,X);
      qq = ((mu_0*Fg_domain)-m)/s;
      v = P.t_distribution_v;
      p = tpdf(qq,v);
    else
      error('unknown P.Fg_distribution_type');
    end
    Fg_prob = p/sum(p);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  else                                                   % GIVEN TRUE Fg
    Fg_domain = M.gene.true_Fg(g);
    nfg = 1;
    Fg_prob = 1;
  end
  M.gene.Fg_post(g,:) = as_row(Fg_prob);
  Fg_prob = as_column(Fg_prob);
  timeused(2) = timeused(2) + toc(ttt);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PROJECTION METHOD
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  pCVmaxs = ones(nfg,1); pCVmids = ones(nfg,1); pCVmins = ones(nfg,1); npat_exps = zeros(nfg,1);
  Fg_prob(Fg_prob<0.01) = 0;  % don't bother considering very unlikely Fg's

  % integrate over Fg distribution
  for fgi=1:nfg, Fg = Fg_domain(fgi);
    if Fg_prob(fgi)==0, continue; end

    %
    %1. determine probability of 0 mutations for each effect
    ttt=tic;
    mu_ssnv = M.pat.mu_c * Fg;
    p0_ssnv = bsxfun(@power,1-mu_ssnv,N_signal_ssnv); %we use a geometric distribution
    p0_ssnv = squeeze(prod(p0_ssnv,2)); %marginalize over categories
    %p0_ssnv = reshape(prod(p0_ssnv,2),M.np,3); %marginalize over categories
    if P.analyze_silent_as_signal
      P0 = p0_ssnv;
    else
      mu_ind = M.pat.mu_ind * Fg;
      p0_ind = (1-mu_ind).^N_signal_indel;
      P0 = [p0_ssnv p0_ind];
    end

    %
    %2. Map effects to degrees for each sample
    P0(isnan(P0))=1; P0(P0>1)=1; P0(P0<0)=0;
    P1 = 1-P0; %probability of seeing at least one mutation in each effect
    [tmp priority] = sort(P1,2,'descend'); shft = (priority - repmat(1:ndeg,M.np,1));
    map = reshape(1:(M.np*ndeg),M.np,ndeg);
    newmap = map + shft*M.np;
    P0 = P0(newmap);
    P1 = P1(newmap);
    has_mutation_sorted = has_mutation(newmap);

    %
    %3. Calculate probability of seeing mutation of degree n
    %   for each sample (+ no mutations at all for degree 0)
    Pdeg = nan(M.np,ndeg+1);
    for d=0:ndeg,
      Pdeg(:,d+1) = prod(P0(:,d+1:end),2);
      if d>0, Pdeg(:,d+1) = Pdeg(:,d+1) .* P1(:,d); end
    end

    %expected number of patients with at least one mutation of any degree 
    npat_exps(fgi) = sum(1-Pdeg(:,1));

    %
    %4. Convert degrees to scores
    Sdeg = [zeros(M.np,1) -log10(P1)];
    Sdeg = min(ceil(Sdeg * 10), 1000);

    %
    %5. Calculate observed degree/score for each patient 
    degree = zeros(M.np,1); score_obs = 0;
    for p = 1:M.np,
      for d = ndeg:-1:1,
	if has_mutation_sorted(p,d), degree(p) = d; break; end
      end
      score_obs = score_obs + Sdeg(p,degree(p)+1);
      M.gene.S(g,p) = Sdeg(p,degree(p)+1);
    end
    score_obs = score_obs / M.gene.min_effect_size(g);
    timeused(3) = timeused(3) + toc(ttt);

    %
    %6. Convolve
    if ~P.skip_pCV_convolutions
      ttt = tic;
      if P.use_pCV_gaussian_approximation && npat_exps(fgi)>=P.pCV_gaussian_approximation_min_npat_exp
        [p m s] = projection_1d_gaussapprox(Sdeg,Pdeg,score_obs);
        M.gene.pCV_Gaussian(g,fgi) = p;
        M.gene.pCV_Gaussian_m(g,fgi) = m;
        M.gene.pCV_Gaussian_s(g,fgi) = s;
        M.gene.pCV_score_obs(g,fgi) = score_obs;
        M.gene.pCV_mefsz(g,fgi) = M.gene.min_effect_size(g);
        pCVmax = p;
        pCVmin = p;
        pCVmid = p;
      else    % convolutions
        numbins = ceil(score_obs+max(5,0.2*score_obs)); H = zeros(numbins,1); newH = zeros(numbins,ndeg+1);
        [pCVmax pCVmin] = projection_1d_convolutions_fast(Sdeg,Pdeg,score_obs,numbins,H,newH);
        pCVmid = pCVmin+rand*(pCVmax-pCVmin);
      end
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % GAUSSIAN BLURRING OF THE NULL DISTRIBUTION
      % --> to be added here
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      pCVmaxs(fgi) = pCVmax; pCVmids(fgi) = pCVmid; pCVmins(fgi) = pCVmin;
      timeused(4) = timeused(4) + toc(ttt);
    end
  end % next member of Fg distribution
  pCVmax = sum(pCVmaxs.*Fg_prob); pCVmid = sum(pCVmids.*Fg_prob); pCVmin = sum(pCVmins.*Fg_prob);
  M.gene.npat_exp(g) = sum(npat_exps.*Fg_prob);

  end  % end of "else" section of "if P.use_v3_test4_method"

  M.gene.pCVmax(g) = pCVmax; M.gene.pCVmid(g) = pCVmid; M.gene.pCVmin(g) = pCVmin;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PREPARATION FOR PERMUTATIONS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  while(true) % "CL+FN container loop", will only be executed once, but we can "break" out at any point to abort going into the permutations

    if P.skip_permutations, break; end

    ttt = tic;

    % find the target regions for this gene
    tidx = find(M.targ.gene_idx==g);
    genelength = sum(M.targ.len(tidx));
    if genelength==0, break; end        % no targets found
  
    % find the mutations for this gene and map to targets
    midx = find(M.mut.gene_idx==g & M.effect.include_in_permutations(M.mut.effect_idx));
    if P.remove_ULTRA_mutations_early
      fprintf('remove_ULTRA_mutations_early\n');
      midx(M.pat.exclude_from_pCL(M.mut.pat_idx(midx)))=[];
    end
    exonstart=1;
    for ti=1:length(tidx),t=tidx(ti);
      z = find(M.mut.chr(midx)==M.targ.chr(t) & M.mut.pos(midx)>=M.targ.start(t) & M.mut.pos(midx)<=M.targ.end(t));
      M.mut.trackpos(midx(z)) = exonstart + M.mut.pos(midx(z))-M.targ.start(t);
      exonstart=exonstart+M.targ.len(t);
    end
    midx(isnan(M.mut.trackpos(midx)))=[];   % remove mutations that didn't map to a target
  
    if P.restrict_to_one_mutation_per_patient   % (chosen randomly)
      midx = midx(randperm(length(midx))); [u ui uj] = unique(M.mut.pat_idx(midx)); midx = midx(ui);
    end
    nm = length(midx);
    M.gene.nmut(g) = nm;
    if nm<2, break; end  % we only do permutations on genes with >=2 mutations

    % read conservation, coverage, and context_and_effect for these regions
    %conservation_track = double(M.FWB.conservation.get(M.targ.chr(tidx),M.targ.start(tidx),M.targ.end(tidx)));
    conservation_track = ones(genelength, 1);
%    conservation_track(conservation_track==200) = NaN;  % missing data
    context_and_effect_track = double(M.FWB.context_and_effect.get(M.targ.chr(tidx),M.targ.start(tidx),M.targ.end(tidx)));
    if isfield(M.FWB,'basewise_coverage')
      coverage_track = double(M.FWB.basewise_coverage.get(M.targ.chr(tidx),M.targ.start(tidx),M.targ.end(tidx)));
    else
      coverage_track = ones(genelength,1);
    end

    % impute at least median coverage in radius around each mutation
    medcov = median(coverage_track(coverage_track>0));
    for i=1:nm,m=midx(i);
      i1 = max(1,M.mut.trackpos(m)-P.radius_to_impute_coverage_around_mutations);
      i2 = min(genelength,M.mut.trackpos(m)+P.radius_to_impute_coverage_around_mutations);
      coverage_track(i1:i2) = max(coverage_track(i1:i2),medcov);
    end
  
    % simplify coverage track (reduce number of discrete coverage levels)
    maxcov = max(coverage_track);
    if maxcov > P.max_coverage_bins
      coverage_track_factor = maxcov / P.max_coverage_bins;
      coverage_track = round(coverage_track / coverage_track_factor); % maybe should be ceil?
      maxcov = max(coverage_track);
    end
    if maxcov==0, break; end   % this gene has no coverage

    % identify subset of the mutations to be included in pCL
    include_in_CL = ~M.pat.exclude_from_pCL(M.mut.pat_idx(midx));

    % enumerate throwable positions for each mutation flavor
    M.effect.permutations_effect_idx(M.effect.permutations_effect_idx < 4) = 1;
    categ = M.mut.categ_idx(midx);
    effect = M.effect.permutations_effect_idx(M.mut.effect_idx(midx));
    categ(effect==4) = inf;  % collapse all indels to categ=inf
    [uce tmp throw_flavor] = unique([categ effect],'rows');  
    nflavors = size(uce,1);
    flavor_counts = histc(throw_flavor,1:nflavors);
    flavor_counts_CL = histc(throw_flavor(include_in_CL),1:nflavors);
    throwable = cell(nflavors,1);
    for fi=1:nflavors
      c = uce(fi,1); e = uce(fi,2);
      if e<4
	%get contexts + effects for this flavor
	ce = find(cellfun(@(x) any(x == c), M.context_and_effect.categ) & ...
	          M.context_and_effect.permutations_effect_idx >= e);
	%ce = find(cellfun(@(x) any(x == c), M.context_and_effect.categ));

	%enumerate positions, weighted by coverage 
	%throwable{fi} = find(ismember(context_and_effect_track, ce));

	throwpos = find(ismember(context_and_effect_track, ce));
	throwable{fi} = cell(length(throwpos), 1);
	for p = [throwpos'; 1:length(throwpos)],
	  throwable{fi}{p(2)} = repmat(p(1), coverage_track(p(1)), 1);
	end
	throwable{fi} = cat(1, throwable{fi}{:});

	%TODO: accommodate allele-specific effects

%        [context_and_effect_idx newbase_idx] = find(M.context_and_effect.Q(:,c,:)==e);
%        cn = [context_and_effect_idx newbase_idx];
%        % make sure the classes of the observed positions are included (in case track has disagreements)
%        mm = midx(throw_flavor==fi);
%        cn_obs = [context_and_effect_track(M.mut.trackpos(mm)) M.mut.newbase_idx(mm)];
%        cn = unique([cn;cn_obs],'rows');
%        % make list of throwable for this flavor
%        nx = size(cn,1); x = cell(nx,1);
%        for i=1:nx
%          idx = find(context_and_effect_track==cn(i,1));
%          y = cell(length(idx),1);
%          for j=1:length(idx), y{j} = repmat([idx(j) cn(i,2)],coverage_track(idx(j)),1); end
%          x{i} = cat(1,y{:});
%        end
%        throwable{fi} = cat(1,x{:});
      else            % e=5=indel: can throw to whole territory, weighted by the coverage track
        x = cell(genelength,1);
	for i=1:genelength, x{i}=i*ones(coverage_track(i),1); end
%        x = cat(1,x{:});
%        throwable{fi} = [x 5*ones(size(x))]; % newbase=5
	throwable{fi} = cat(1, x{:});
      end
    end

    nthrowable = cellfun(@(x) size(x,1),throwable);
    if sum(nthrowable)==0, break; end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ANALYTICAL CLUSTERING METHOD by Julian Hess
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % pCL2
    x0 = max(histc(M.mut.trackpos(midx(include_in_CL)), 1:genelength));
    pclust2 = max(1e-16,pCLmn(x0, nm, genelength));

    % pFN2 and pFN3 (Gaussian approximation)
    fcidx = find(nthrowable > 0)'; %we only convolve in throwable areas.
    m=[]; s=[];
    for i = fcidx
      cons = conservation_track(throwable{i}(:, 1));
      n = flavor_counts(i);
      m(end+1,1) = n*mean(cons); s(end+1,1) = sqrt(n*max(0.1,std(cons)).^2);
      tmp = histc(cons, 0:100); tmp = tmp./sum(tmp);
      tmp2 = iconv(tmp, flavor_counts(i));
      if i == fcidx(1), cons_PDF_conv = tmp2;
      else cons_PDF_conv = conv(cons_PDF_conv, tmp2); end
    end
    m = sum(m); s = sqrt(sum(s.^2));
    ocons = nansum(conservation_track(M.mut.trackpos(midx)));
    if isnan(ocons), ocons = min_cons*sum(flavor_counts); end
    if ocons>length(cons_PDF_conv)-1, pcons2 = 1; else pcons2 = max(1e-16,1-sum(cons_PDF_conv(1:ocons))); end
    M.gene.pFN_Gaussian_m(g) = m;
    M.gene.pFN_Gaussian_s(g) = s;
    M.gene.pFN_ocons(g) = ocons;

    % pCLFN2
    pjoint2 = fisher_combined_p([pclust2 pcons2]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ORIGINAL PERMUTATIONS METHOD by Petar Stojanov
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % range of metrics and bins for joint distrib
    min_clust = 0; max_clust = 1;
    %min_cons = min(conservation_track); max_cons = max(conservation_track);
    %if isnan(min_cons), min_cons = 0; end
    %if isnan(max_cons), max_cons = 100; end
    nbins = 100;
    binsize_clust = (max_clust-min_clust)/(nbins-1);
    %binsize_cons = (max_cons-min_cons)/(nbins-1);
    
    % calculate metrics for the observed mutations
    %obs_cons = nanmean(conservation_track(M.mut.trackpos(midx)));
    %if isnan(obs_cons), obs_cons = min_cons; end
    %obs_cons = obs_cons / P.pFN_min_effect_size;
    %obs_cons_bin = 1+floor((obs_cons - min_cons)/binsize_cons);
    obs_clust = new_clustering_statistic(M.mut.trackpos(midx(include_in_CL)),genelength,P.clustering_metric,true) / P.pCL_min_effect_size;
    if isnan(obs_clust), obs_clust = min_clust; end
    obs_clust_bin = 1+floor((obs_clust - min_clust)/binsize_clust);

    timeused(5) = timeused(5) + toc(ttt);
  
    %%%%%%%%%%%%%%%%%%%
    % PERMUTATIONS
    %%%%%%%%%%%%%%%%%%%

    % make mask for permuted position inclusion in pCL
    perm_include_in_pCL = cell(nflavors,1);
    for j=1:nflavors, if nthrowable(j)>0, perm_include_in_pCL{j} = [true(flavor_counts_CL(j),1);false(flavor_counts(j)-flavor_counts_CL(j),1)]; end,end
    perm_include_in_pCL = cat(1,perm_include_in_pCL{:});
    all_excluded_in_pCL = all(~perm_include_in_pCL);

    ttt = tic;
    k_clust = 0; %k_cons = 0; joint_hist = zeros(nbins,nbins);
    nperm = 0; first_check = 100; check_every = 10000; finished = false; tt1=tic;
    while(~finished), nperm = nperm + 1;
      % randomly throw mutations
      thrown_muts = cell(nflavors,1);
      for j=1:nflavors, if nthrowable(j)>0, thrown_muts{j} = throwable{j}(ceil(nthrowable(j)*rand(flavor_counts(j),1)),:); end,end
      perm_mutpos = cat(1,thrown_muts{:});

      % calculate metrics for this permutation and increment histograms
      %perm_cons = nanmean(conservation_track(perm_mutpos(:,1)));
      %if isnan(perm_cons), perm_cons = max_cons; end
      perm_clust = new_clustering_statistic(perm_mutpos(perm_include_in_pCL,1),genelength,P.clustering_metric);
      if isnan(perm_clust) || all_excluded_in_pCL, perm_clust = max_clust; end
      %if perm_cons>=obs_cons, k_cons=k_cons+1; end
      if perm_clust>=obs_clust, k_clust=k_clust+1; end
      %bin_cons = 1+floor((perm_cons - min_cons)/binsize_cons);
      bin_clust = 1+floor((perm_clust - min_clust)/binsize_clust);
      %joint_hist(bin_clust,bin_cons) = joint_hist(bin_clust,bin_cons) + 1;
      
      if nperm~=first_check && mod(nperm,check_every)>0 && nperm<P.maxperm, continue; end  % don't need to check p-values every iteration
        
      % calculate marginal and joint p-values
      [pclust ci_ratio_clust] = calc_pval_and_ci_ratio(k_clust,nperm);
%      [pcons ci_ratio_cons] = calc_pval_and_ci_ratio(k_cons,nperm);
%      landfilled_hist = landfill(joint_hist/nperm);
%      bin_score = -log10(landfilled_hist);
%      obs_score = bin_score(obs_clust_bin,obs_cons_bin);
%      k_joint = sum(joint_hist(bin_score>=obs_score));
%      [pjoint ci_ratio_joint] = calc_pval_and_ci_ratio(k_joint,nperm);
%      max_ci_ratio = max([ci_ratio_joint,ci_ratio_cons,ci_ratio_clust]);
      finished = (ci_ratio_clust<=P.theta) | (nperm>=P.maxperm);

      % cap p-values
      pcap = 10^-nm;  %pcons = max(pcons,pcap);
      pclust = max(pclust,pcap); %pjoint = max(pjoint,pcap);
      pcap = 1/nperm; %pcons = max(pcons,pcap);
      pclust = max(pclust,pcap); %pjoint = max(pjoint,pcap);
      
      % STATUS REPORT
      sec_remaining = toc(tt) * (length(genes_to_analyze)-num_genes_analyzed)/max(1,num_genes_analyzed);
      eta = now + (sec_remaining/(60*60*24));
      timeprofile = sprintf('%.0f ',100*timeused/sum(timeused));
      tmp = {'%s %s %5d/%5d %-15s %-4d %7d  CV %-0.16f CL %-0.16f FN %-0.16f  [ %-4d %-4d %-4d %-4d %-4d %-4d ]  %s\n',...
             datestr(now,'dddHH:MM:SS'),datestr(eta,'dddHH:MM:SS'),g,M.ng,M.gene.name{g},nm,nperm,pCVmax,pclust,pcons,...
             qqhist(1),qqhist(2),qqhist(3),qqhist(4),qqhist(5),qqhist(6),timeprofile};
      fprintf(tmp{:});
      gene_status_reported = true;
      
    end  % next permutation

    % record results of permutations
    M.gene.sCL(g) = obs_clust;
    %M.gene.sFN(g) = obs_cons;
    M.gene.nperm(g) = nperm;
    M.gene.pCL(g) = pclust;
    %M.gene.pFN(g) = pcons;
    %M.gene.pCLFN(g) = pjoint;
    M.gene.pCL2(g) = pclust2;
    M.gene.pFN2(g) = pcons2;
    M.gene.pCLFN2(g) = pjoint2;
    timeused(6) = timeused(6) + toc(ttt);
    break
  end   % "CL+FN container loop"

  num_genes_analyzed = num_genes_analyzed + 1;

  if ~gene_status_reported
      % STATUS REPORT
      sec_remaining = toc(tt) * (length(genes_to_analyze)-num_genes_analyzed)/max(1,num_genes_analyzed);
      eta = now + (sec_remaining/(60*60*24));
      timeprofile = sprintf('%.0f ',100*timeused/sum(timeused));
      tmp = {'%s %s %5d/%5d %-15s %-4d %7d  CV %-0.16f CL %-0.16f FN %-0.16f  [ %-4d %-4d %-4d %-4d %-4d %-4d ]  %s\n',...
             datestr(now,'dddHH:MM:SS'),datestr(eta,'dddHH:MM:SS'),g,M.ng,M.gene.name{g},nm,nperm,pCVmax,pclust,pcons,...
             qqhist(1),qqhist(2),qqhist(3),qqhist(4),qqhist(5),qqhist(6),timeprofile};
      fprintf(tmp{:});
      gene_status_reported = true;
  end

  % increment text representation of q-q plot for status report
  sp = min(6,floor(-log10(min([pCVmax pclust pcons]))));
  if sp>=1, qqhist(sp)=qqhist(sp)+1; end

  catch me
   fprintf('ERROR with gene %s\n',M.gene.name{g});
   disp(me.message)
   disp(me)
   for i=1:length(me.stack), disp(me.stack(i)); end

  end

end, fprintf('\n');   % next gene

% close track files
if ~P.skip_permutations,
  M.FWB.context_and_effect.close();
end

% finish pFN3 calculation
M.gene.pFN3 = 1-normcdf(M.gene.pFN_ocons,M.gene.pFN_Gaussian_m,M.gene.pFN_Gaussian_s);

% get rid of NaNs in pCV
M.gene.pCVmin(isnan(M.gene.pCVmin)) = 1;
M.gene.pCVmid(isnan(M.gene.pCVmid)) = 1;
M.gene.pCVmax(isnan(M.gene.pCVmax)) = 1;

% cap p-values
pcap = 1e-16; M.gene.pCVmin = max(M.gene.pCVmin,pcap); M.gene.pCVmid = max(M.gene.pCVmid,pcap); M.gene.pCVmax = max(M.gene.pCVmax,pcap);

% combine p-values
M.gene.pmid = max(1e-16,fisher_combined_p([M.gene.pCVmid M.gene.pCLFN]));
M.gene.pmax = max(1e-16,fisher_combined_p([M.gene.pCVmax M.gene.pCLFN]));

% for genes that didn't undergo permutations, just take the MutSigCV values
%idx = find(isnan(M.gene.pCLFN)|isnan(M.gene.nperm)|M.gene.nperm==0);
%M.gene.nperm(idx)=0;
%M.gene.pFN(idx)=nan; M.gene.pCL(idx)=nan; M.gene.pCLFN(idx)=nan;
%M.gene.pmid(idx) = M.gene.pCVmid(idx); M.gene.pmax(idx) = M.gene.pCVmax(idx);

% FDR
M.gene.qCF = calc_fdr_value(M.gene.pCF);
M.gene.qCV = calc_fdr_value(M.gene.pCVmax);
M.gene.qCL = calc_fdr_value(M.gene.pCL);
M.gene.qFN = calc_fdr_value(M.gene.pFN);
M.gene.qCL2 = calc_fdr_value(M.gene.pCL2);
M.gene.qFN2 = calc_fdr_value(M.gene.pFN2);
M.gene.q = calc_fdr_value(M.gene.pmax);

timeprofile = sprintf('%.0f ',100*timeused/sum(timeused));
fprintf('Finished in %f sec [%s]\n',toc(tt),timeprofile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% RESULTS 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = M.gene;
G = sort_struct(G,{'pmax','npat'},[1 -1]);
%G = rename_field(G,{'name','pCVmax','pmax','nsyn','nnon','Nsyn','Nnon'},{'gene','pCV','p','nsil','nstp','Nsil','Nstp'});
%G.nnon = G.nmis+G.nstp+G.nspl+G.nind;
G.nnon = sum([G.nssnv(:, ~M.effect.is_background & ~M.effect.is_indel) G.nind], 2);
G.freq = G.npat/M.np;
G.freq_excess = max(0,(G.npat-G.npat_exp))/M.np;
G = move_field_to_after(G,'freq','nsite');
G = move_field_to_after(G,'freq_excess','freq');

if P.skip_output_files, return; end

nsg_CF = sum(G.qCF<=0.1);
nsg_CV = sum(G.qCV<=0.1);
nsg_CL = sum(G.qCL<=0.1);
nsg_FN = sum(G.qFN<=0.1);
nsg_CL2 = sum(G.qCL2<=0.1);
nsg_FN2 = sum(G.qFN2<=0.1);
nsg = sum(G.q<=0.1);
fprintf('genes with q<=0.1:   pCV %-5d pCL %-5d pFN %-5d pCL2 %-5d pFN2 %-5d pCF %-5d total %-5d\n',nsg_CV,nsg_CL,nsg_FN,nsg_CL2,nsg_FN2,nsg_CF,nsg);

%memaudit

fprintf('Saving results... ');

% save patient counts and rates file
pat=M.pat;
save_textfile(num2str(M.np),[outdir '/num_patients.txt']);
pat = rmfield_if_exist(pat,{'N_c','n_c','N_c_2','mu_c'});
save_struct(pat,[outdir '/patient_counts_and_rates.txt']);
save([outdir '/patients.mat'],'pat');

% save categories file
save_struct(M.categ,[outdir '/mutcategs.txt']);

% save full output matfile
if ~P.skip_verbose_output_mat
  save([outdir '/results.verbose.PARTIAL.mat'],'G');
  movefile([outdir '/results.verbose.PARTIAL.mat'],[outdir '/results.verbose.mat']);
end
% and version with large fields removed
G = rmfield_if_exist(G,{'bagel','S'});
save([outdir '/results.PARTIAL.mat'],'G');
movefile([outdir '/results.PARTIAL.mat'],[outdir '/results.mat']);

% save sig_genes table
g = keep_fields_if_exist(G,{'gene','longname','codelen','nnei','nncd','nsil','nmis','nstp','nspl','nind',...
                    'nnon','npat','nsite','pCV','pCL','pFN','pCL2','pFN2','pCF','p','q'});
save_struct(g,[outdir '/sig_genes.PARTIAL.txt']);
movefile([outdir '/sig_genes.PARTIAL.txt'],[outdir '/sig_genes.txt']);

% number of significant genes
save_textfile(num2str(nsg),   [outdir '/num_sig_genes.txt']);
save_textfile(num2str(nsg_CV),[outdir '/num_sig_genes.pCV.txt']);
save_textfile(num2str(nsg_CL),[outdir '/num_sig_genes.pCL.txt']);
save_textfile(num2str(nsg_FN),[outdir '/num_sig_genes.pFN.txt']);
save_textfile(num2str(nsg_CL2),[outdir '/num_sig_genes.pCL2.txt']);
save_textfile(num2str(nsg_FN2),[outdir '/num_sig_genes.pFN2.txt']);
save_textfile(num2str(nsg_CF),[outdir '/num_sig_genes.pCF.txt']);

% save final mutation list (do last because it takes the most time)
if ~P.skip_output_maf
  fprintf('Saving final_analysis_set.maf...\n');
  M.mut.categ = M.mut.categ_idx;
  save_struct(M.mut,[outdir '/final_analysis_set.PARTIAL.maf']);
  movefile([outdir '/final_analysis_set.PARTIAL.maf'],[outdir '/final_analysis_set.maf']);
  save2(M.mut,[outdir '/final_analysis_set.PARTIAL.mafdir']);
  movefile([outdir '/final_analysis_set.PARTIAL.mafdir'],[outdir '/final_analysis_set.mafdir']);
end

save_textfile('done',[outdir '/finished.txt']);
fprintf('Done.  Results written to %s\n',outdir);

% close logfile and msgfile
fclose(msgfile);
fclose(logfile);

%memaudit


%  function memaudit
%    if P.memory_audit
%      fprintf('MEMAUDIT @ ');
%      disp(datestr(now));
%      whos
%    end
%  end
%
%  function logprintf(varargin)
%    txt = sprintf(varargin{:});
%    fprintf('%s',txt);
%    fprintf(logfile,'%s',txt);
%  end
%
%  function fprintf(varargin)
%    txt = sprintf(varargin{:});
%    fprintf('%s',txt);
%    fprintf(msgfile,'%s',txt);
%  end
    

end





