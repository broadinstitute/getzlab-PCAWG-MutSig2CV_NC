function [M,P] = MutSig2CV_v5_load(mutation_file, outdir, effect_table, coverage_models_mat_file, ...
                                         context_and_effect_fwb_file, context_and_effect_categs_file, ...
                                         covariates_file, FixedWidthBinary_jar_file, params_file)

if nargin<9, error('requires 9 input arguments'); end
if nargin>10, error('extraneous arguments'); end

% ensure input files exist
if ischar(mutation_file), demand_file(mutation_file); end
if ~strcmp('',effect_table), demand_file(effect_table); end
demand_file(coverage_models_mat_file);
demand_file(context_and_effect_fwb_file);
demand_file(context_and_effect_categs_file);
demand_file(covariates_file);
demand_file(FixedWidthBinary_jar_file);

% add jar to java classpath
javaclasspath(FixedWidthBinary_jar_file);

% ensure output directory OK
if ~exist('outdir','var') || isempty(outdir) || strcmpi('none',outdir), outdir = ['/tmp/mutsig_temp' num2str(round(now*1e6+1000*rand))]; %fprintf('Writing to temporary outdir %s\n',outdir);
end
if ~isempty(outdir)
  ede(outdir);
  if exist([outdir '/results.mat'],'file')
    fprintf('WARNING: results already exist in this directory!\n');
    fprintf('Pausing 10 seconds before continuing.\n');
    pause(10);
  end
end
%msgfile = fopen([outdir '/LOAD_msg.txt'],'wt');

P=[];
if exist('params_file','var') && ~strcmp(params_file,'') && ~strcmp(params_file,'none')
  P = process_params_file(P,params_file);
end

% MutSig LOAD parameters
P = impose_default_value(P,'scan_for_and_remove_duplicate_patients',true);
P = impose_default_value(P,'enforce_target_list',true);
P = impose_default_value(P,'enforce_target_list_inclusive_of_noncoding',true);
P = impose_default_value(P,'enforce_effect',true);
P = impose_default_value(P,'use_new_find_mut_categs',true);
P = impose_default_value(P,'number_of_categories_to_discover',8);
P = impose_default_value(P,'categs_file','');
if P.number_of_categories_to_discover==96
  if ~isempty(P.categs_file), error('please don''t specify a categs_file if using 96 categories'); end
end
if isfield(P,'callschemes_all_WGS'), error('please use P.force_all_wgs_callschemes'); end
P = impose_default_value(P,'force_all_wgs_callschemes',false);
if P.force_all_wgs_callschemes, fprintf('Will assume all WGS callschemes.\n'); end
P = impose_default_value(P,'callschemes_file','');
P = impose_default_value(P,'complete_missing_fields',false);
P = impose_default_value(P,'strand_specific',false);
P = impose_default_value(P,'gene_min_frac_coverage_required',0.10);
P = impose_default_value(P,'memory_audit',false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LOAD   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%memaudit

fprintf('LOADING DATA\n');

loadstart = tic();

M=[];
if ischar(mutation_file), M.maf = mutation_file; end

% FWB tracks
M.FWB = [];
M.FWB.context_and_effect_file = context_and_effect_fwb_file;
demand_file(context_and_effect_fwb_file);
M.context_and_effect = load_struct(context_and_effect_categs_file);
M.context_and_effect.context65 = map_categories_to_65(context_and_effect_categs_file); %XXX: make sure this works

% open FWB tracks
M.FWB.context_and_effect = org.broadinstitute.cga.tools.seq.FixedWidthBinary(M.FWB.context_and_effect_file);

%memaudit

% MUTATIONS 
fprintf('Loading mutations   ');
if isstruct(mutation_file)
  fprintf('(from passed struct)  ');
  M.mut = mutation_file;
elseif exist(mutation_file,'dir')
  % if input file ends in .M or .Mdir just try to loadM it and then we're done!
  if grepm('\.M$',{mutation_file}) | grepm('\.Mdir$',{mutation_file})
    fprintf('(from saveM)  ');
    tic;M = loadM(mutation_file);toc
    demand_fields(M,{'mut','pat','gene','cov','effect'});
    return
  else
    % it's just a mutation save2 file, try to load2 it
    fprintf('(from save2)  ');
    tic; M.mut = load2(mutation_file);toc
  end
elseif exist(mutation_file,'file')
  fprintf('(from file)  ');
  tic;M.mut = load_struct(mutation_file);toc
end

%memaudit

% compute newbase if necessary
if isfield(M.mut,'newbase') && ~any(cellfun('isempty',M.mut.newbase))
  % already ok
elseif isfield(M.mut,'Tumor_Seq_Allele1')
  M.mut.newbase = M.mut.Tumor_Seq_Allele1;
  if isfield(M.mut,'Reference_Allele') && isfield(M.mut,'Tumor_Seq_Allele2')
    idx = find(strcmp(M.mut.Reference_Allele,M.mut.Tumor_Seq_Allele1));
    M.mut.newbase(idx) = M.mut.Tumor_Seq_Allele2(idx);
  end
elseif isfield(M.mut,'tum_allele1')
  M.mut.newbase = M.mut.tum_allele1;
  if isfield(M.mut,'ref_allele') && isfield(M.mut,'tum_allele2')
    idx = find(strcmp(M.mut.ref_allele,M.mut.tum_allele1));
    M.mut.newbase(idx) = M.mut.tum_allele2(idx);
  end
elseif isfield(M.mut,'Tumor_Seq_Allele2')
  M.mut.newbase = M.mut.Tumor_Seq_Allele2;
elseif isfield(M.mut,'tum_allele2')
  M.mut.newbase = M.mut.tum_allele2;
end

% if P.complete_missing_fields, M.mut = complete_maf(M.mut); end

% make sure required fields exist (renaming if necessary)
rf = {
    {'patient','Tumor_Sample_Barcode','Patient_name'},
    {'chr','Chromosome'},
    {'pos','Position','start','Start_position'},
    {'ref_allele','Reference_Allele'},
    {'newbase','Tumor_Allele','Tum_allele','Alt_allele','Alternate_allele','Tumor_Seq_Allele2'},
};
if ~P.enforce_target_list
  rf = [rf;{
    {'gene','Hugo_Symbol','Gene_name'},
  }];
end
if ~P.enforce_effect
  rf = [rf;{
    {'type','Variant_Classification'};
    {'classification','Variant_Type'};
    }];
end

f = fieldnames(M.mut);
for i=1:length(rf)
  matches = find(ismember(lower(f),lower(rf{i})));
  if isempty(matches)
    fprintf('\nMutation file is missing a column named one of the following:\n');
    pr(rf{i});
    error('Mutation file missing %s column.',rf{i}{1});
  end
  if length(matches)>1
    fprintf('\nMutation file contains multiple columns for %s info:\n', rf{i}{1});
    pr(f(matches));
    match_to_first = find(strcmpi(lower(f),lower(rf{i}{1})),1);
    if ~isempty(match_to_first), matches = match_to_first; else matches = matches(1); end
    fprintf('Will use %s\n',f{matches});
  end
  M.mut = rename_field(M.mut,f{matches},rf{i}{1});
end

%memaudit

% convert chr/pos
fprintf('Parsing mutation positions...\n');
M.mut.chr = convert_chr(M.mut.chr);
h = histc(M.mut.chr,1:24);
if all(h(1:5)>=10) && all(h(1:22)>=1) && h(23)==0
  fprintf('\n');
  fprintf('\t***********************************************\n');
  fprintf('\t*                                             *\n');
  fprintf('\t*                 W A R N I N G               *\n');
  fprintf('\t*                                             *\n');
  fprintf('\t*   Looks like chrX mutations got removed!    *\n');
  fprintf('\t*                                             *\n');
  fprintf('\t* Possible cause:                             *\n');
  fprintf('\t* make_numeric(chr) -> str2double(''X'') -> NaN *\n');
  fprintf('\t*                                             *\n');
  fprintf('\t***********************************************\n');
  fprintf('\n');
  fprintf('Pausing for 10 minutes.\n');
  pause(600);
end
idx = find(isnan(M.mut.chr));
if ~isempty(idx)
  fprintf('Omitting %d/%d mutations on nonstandard chromosomes\n',length(idx),slength(M.mut));
  M.mut = reorder_struct_exclude(M.mut,idx);
end
M.mut = make_numeric(M.mut,'pos');
idx = find(isnan(M.mut.pos));
if ~isempty(idx)
  fprintf('Omitting %d/%d mutations at indecipherable positions\n',length(idx),slength(M.mut));
  M.mut = reorder_struct_exclude(M.mut,idx);
end
if isfield(M.mut,'i_tumor_f'), M.mut = make_numeric(M.mut,'i_tumor_f'); end

% remove duplicate mutations
tmp=[];
[tmp.u tmp.ui tmp.uj] = unique_combos(M.mut.patient,M.mut.chr,M.mut.pos);
if length(tmp.ui)<slength(M.mut)
  fprintf('Keeping %d/%d unique mutations.\n',length(tmp.ui),slength(M.mut));
  M.mut = reorder_struct(M.mut,tmp.ui);
end

%memaudit

if P.scan_for_and_remove_duplicate_patients
  % remove duplicate patients
  fprintf('Scanning for duplicate patients...\n');
  tmp = new_find_duplicate_samples(M.mut,P);
  if ~isempty(tmp.drop)
    fprintf('Removing the following %d duplicate patients:\n',length(tmp.drop));
    disp(tmp.drop);
    M.mut = reorder_struct_exclude(M.mut,ismember(M.mut.patient,tmp.drop));
  end
  tmp = [];
end

%memaudit

% COVERAGE MODELS
fprintf('Loading coverage models...\n');
tmp = load(coverage_models_mat_file,'C');
demand_fields(tmp.C.cov,{'CE_cov','CE_terr'});
if any(isinf(tmp.C.cov.CE_cov(:))) || any(isinf(tmp.C.cov.CE_terr(:))), error('infs in terr/cov'); end
if any(isnan(tmp.C.cov.CE_cov(:))) || any(isnan(tmp.C.cov.CE_terr(:))), error('nans in terr/cov'); end
M.cov = tmp.C.cov;
M.gene = tmp.C.gene;

if isfield(tmp.C, 'targ'),
  M.targ = tmp.C.targ;
end

tmp=[];

%memaudit

% TARGET LIST
fprintf('Processing target list.\n');
if isfield(M, 'targ'),
  fprintf('Using list embedded in coverage model.\n');
else
  error('Malformed coverage model -- no target list present!\n');
end
demand_fields(M.targ,{'gene','chr','start','end'});
M.targ.chr = convert_chr(M.targ.chr);
M.targ = make_numeric(M.targ,{'start','end'});
M.targ = reorder_struct_exclude(M.targ,isnan(M.targ.chr)|isnan(M.targ.start)|isnan(M.targ.end));
if any(M.targ.start>M.targ.end), error('target list has start>end'); end
M.targ = sort_struct(M.targ,{'chr','start','end'});
M.targ.len = M.targ.end-M.targ.start+1;
%if mean(M.targ.len)>5000, error('Looks like target list is whole transcripts.  Need individual exons!'); end

% enforce target list (if requested)
%XXX: use mm2t fast?  if we want to consider flanks, use iterative procedure as performed in India prep
if P.enforce_target_list
  fprintf('Enforcing target list.\n');
  %PP=P;
  %PP.enforce_target_list_footprint_table = M.cov.gene; % (only relevant if P.enforce_target_list_inclusive_of_noncoding==true)
  tmp=[]; %tmp.map = map_mutations_to_targets(M.mut,M.targ,PP);

  M.mut = sort_struct(M.mut, {'chr' 'pos'});
  tmp.map = map_mutations_to_targets_fast([M.mut.chr M.mut.pos], [M.targ.chr M.targ.start M.targ.end]);

  idx = find(tmp.map == -1);
  if ~isempty(idx)
    fprintf('Removing %d/%d mutations that fall outside target gene intervals.\n',length(idx),slength(M.mut));
    M.mut = reorder_struct_exclude(M.mut,idx);
    tmp.map(idx)=[];
  end
  tmp.new_gene = M.targ.gene(tmp.map);
  if isfield(M.mut,'gene')
    tmp.old_gene = M.mut.gene;
    idx = find(~strcmpi(tmp.old_gene,tmp.new_gene));
    if ~isempty(idx)
      fprintf('Reassigning the following %d gene identities:', length(idx));
      tmp.reassign = stringsplice([tmp.old_gene(idx) tmp.new_gene(idx)],1,' -> ');
      count(tmp.reassign,1);
    end
  end
  M.mut.gene = tmp.new_gene;
end
tmp=[];

% remove IGR mutations
idx = find(strcmp('',M.mut.gene)|strcmpi('unknown',M.mut.gene));
if ~isempty(idx)
  fprintf('Omitting %d/%d mutations that are deep-intergenic/unknown-gene.\n',length(idx),slength(M.mut));
  M.mut = reorder_struct_exclude(M.mut,idx);
end

% EFFECT 
% supports generic effects 
if strcmp(effect_table, ''),
  M.effect = [];
  M.effect.name = {'ncd';'syn';'mis';'non';'spl';'indel_ncd';'indel_cod';'indel_spl'};
  M.effect.permutations_effect_idx = [1;2;3;4;4;1;5;5];
  M.effect.include_in_permutations = (M.effect.permutations_effect_idx>=3);
  M.effect.is_indel = grepmi('indel',M.effect.name);
  M.effect.newbase_dependent = [0; 1; 1; 1; 0; 0; 0; 0];
  M.effect.is_background = [1; 1; 0; 0; 0; 1; 0; 0];
else
  M.effect = load_struct(effect_table);
  demand_fields(M.effect, {'name' 'permutations_effect_idx' 'include_in_permutations' 'is_indel' 'newbase_dependent' 'is_background'});
  M.effect.permutations_effect_idx = str2double(M.effect.permutations_effect_idx);
  M.effect = make_logical(M.effect, {'include_in_permutations' 'is_indel' 'newbase_dependent' 'is_background'});
end
M.effect.deg_ord = NaN(slength(M.effect), 1);
M.effect.deg_ord(find(~M.effect.is_background)) = 1:nnz(~M.effect.is_background);

if P.enforce_effect
  fprintf('Enforcing "effect" based on reference files\n');

  %TODO: support enforcing newbase-dependent effects
  %for now, we hardcode something to munge the busted e15 tracks
  M.mut.effect_idx = NaN(slength(M.mut), 1);
 
  %%explicitly annotate indels
  nb_l = cellfun(@length, M.mut.newbase);
  ra_l = cellfun(@length, M.mut.ref_allele);
  M.mut.effect_idx(nb_l > 1 & nb_l > ra_l) = -1;
  M.mut.effect_idx(ra_l > 1 & ra_l > nb_l) = -1;

  %single base indels
  M.mut.effect_idx((strcmp(M.mut.ref_allele, '-') & nb_l == 1) | (strcmp(M.mut.newbase, '-') & ra_l == 1)) = -1;

  M.mut.is_indel = M.mut.effect_idx == -1;

  %%explicitly annotate [DTO]NVs
  onv_idx = nb_l > 1 & nb_l == ra_l;
  M.mut.effect_idx(onv_idx) = -2;

  %
  %map to effect_idx 
  M.mut.context_and_effect = double(M.FWB.context_and_effect.get(M.mut.chr, M.mut.pos));
  M.mut.context_and_effect(M.mut.context_and_effect == -1) = NaN;
  %get context65 while we're at it
  M.mut.context65 = nansub(M.context_and_effect.context65, M.mut.context_and_effect);

  %%sSNVs
  idx = isnan(M.mut.effect_idx);
  M.mut.effect_idx(idx) = mod(M.mut.context_and_effect(idx), 5); M.mut.effect_idx(idx & M.mut.effect_idx == 0) = 5;

  %%ONVs -> take mean and round away from zero
  onv_CE = cell(nnz(onv_idx), 1);
  for x = [find(onv_idx)'; 1:length(onv_CE)],
    i = x(1); j = x(2);
    onv_CE{j} = mod(double(M.FWB.context_and_effect.get(M.mut.chr(i), M.mut.pos(i), ...
                                                                  M.mut.pos(i) + nb_l(i) - 1)), 5);
    onv_CE{j}(onv_CE{j} == 0) = 5;
  end
  %M.mut.effect_idx(onv_idx) = cellfun(@(x) round(mean(ceil(abs(x)).*sign(x - 2))), onv_CE) + 2;
  M.mut.effect_idx(onv_idx) = cellfun(@(x) round(mean(x)), onv_CE);

  %%indel -> get their own effect
  %NOTE: we may want to change this ...
  M.mut.effect_idx(M.mut.is_indel) = 6; %TODO: this shouldn't be hardcoded.

  %XXX: how to handle flanking?
else
  error('Not yet supported.')
%  fprintf('Looking up "effect" in mutation_type_dictionary_file\n');
%  D = load_struct(mutation_type_dictionary_file);
%  demand_fields(D,'effect');
%  f = setdiff(fieldnames(D),{'effect'});
%  demand_fields(M.mut,f);
%  D.effect_idx = listmap(D.effect,M.effect.name);
%  tmp=[]; tmp.di = multimap(M.mut,D,f);
%  M.mut.effect_idx = nansub(D.effect_idx,tmp.di);
end

idx = find(isnan(M.mut.effect_idx));
if ~isempty(idx)
  fprintf('Omitting %d/%d mutations of unknown "effect"\n',length(idx),slength(M.mut));
  pr(M.mut,{'gene','patient','classification','type','chr','pos','ref_allele','newbase'},idx(1:min(1e3,length(idx))))
  if length(idx)>1e3, fprintf('...\n'); end
  M.mut = reorder_struct_exclude(M.mut,idx);
end

% more mutation conversions
fprintf('Annotating mutation effects...\n');
M.mut.newbase_idx = listmap(regexprep(M.mut.newbase,'^(.).*$','$1'), {'A', 'C', 'G', 'T'});
M.mut.newbase_idx(M.effect.is_indel(M.mut.effect_idx)) = 5;  % indels get newbase_idx=5
idx = find(isnan(M.mut.newbase_idx));
if ~isempty(idx)
  fprintf('Omitting %d/%d mutations of invalid "newbase"\n',length(idx),slength(M.mut));
  pr(M.mut,{'gene','patient','classification','type','chr','pos','ref_allele','newbase'},idx(1:min(1e3,length(idx))))
  if length(idx)>1e3, fprintf('...\n'); end
  M.mut = reorder_struct_exclude(M.mut,idx);
end
if isfield(M.mut,'invisible_to_CLFN')
  fprintf('WARNING:  "invisible_to_CLFN"   no longer supported!  This field will be ignored!  Please use   P.coding_BMR_threshold_for_pCL_exclusion   instead.\n');
  M.mut = rmfield(M.mut,'invisible_to_CLFN');
end

% PATIENTS 

[M.pat.name tmp M.mut.pat_idx] = unique(M.mut.patient);
M.np = slength(M.pat);
fprintf('%d patients\n',M.np);
if M.np<2, fprintf('WARNING:  MutSig is not applicable to single patients.\n'); end

%memaudit

% CALLSCHEMES AND TERR/COV
%if ~strcmp(P.callschemes_file, '') & demand_file(P.callschemes_file),
%  % LOAD FROM CALLSCHEMES FILE
%  Cs = load_struct(P.callschemes_file);
%  demand_fields(Cs, {'pat' 'callscheme'}); Cs = make_numeric(Cs, 'callscheme');
%  if any(isnan(listmap(M.pat.name, Cs.pat))), error('Incomplete callscheme specification.'); end
%  if any(~ismember(Cs.callscheme, 0:3)), error('Callschemes must be specified as 0-3.  See impute_callschemes.m for details.'); end
%  callscheme_names = {'coding only';'exome+100bp flanks';'all capture (no interval list)';'all genome (WGS)'};
%  M.pat.callscheme = mapacross(M.pat.name, Cs.pat, Cs.callscheme);
%  M.pat.callscheme_name = callscheme_names(M.pat.callscheme + 1);
%  %also calculate patient counts for counts/rates file
%  M.mut = add_and_convert_simple_fieldnames(M.mut);
%  M.mut = add_helper_is_fields(M.mut); 
%  pvec = 1:slength(M.pat); 
%  M.pat.nmut = as_column(histc(M.mut.pat_idx,1:slength(M.pat)));
%  M.pat.n_coding = as_column(histc(M.mut.pat_idx(M.mut.is_coding), pvec));
%  M.pat.n_flank = as_column(histc(M.mut.pat_idx(M.mut.is_flank), pvec));
%  M.pat.n_coding_nonsilent = as_column(histc(M.mut.pat_idx(M.mut.is_coding & ~M.mut.is_silent), pvec));
%  M.pat.fracflank = M.pat.n_flank./M.pat.nmut;
%else % IMPUTE CALLSCHEMES
%  M = impute_callschemes(M,P); % (consults P.force_all_wgs_callschemes)
%  if ~all(M.pat.callscheme>=0 & M.pat.callscheme<=3), error('Problem imputing callschemes'); end
%end
%
%M.pat.cov_idx = mapacross(M.pat.callscheme,[0 1 2 3],[1 2 3 5]);
%M.pat = rmfield(M.pat,'callscheme');
%count(M.pat.callscheme_name);

% remove genes with very low coverage
%NOTE: this really should be moved up -- should reorder the target list before MM2T.
gene_removal_idx = sum(M.cov.CE_cov, 2)./sum(M.cov.CE_terr, 2) < P.gene_min_frac_coverage_required;

if nnz(idx)
  fprintf('Omitting %d/%d genes because they have extremely low coverage.\n',nnz(idx),size(M.cov.CE_terr, 2));
  M.cov = reorder_struct(M.cov, ~gene_removal_idx);
  M.gene = reorder_struct(M.gene, ~gene_removal_idx);

  %TODO: for posterity, need to reorder the target list too.
end

%memaudit

% GENES
M.ng = slength(M.gene);
M.targ.gene_idx = listmap(M.targ.gene,M.gene.name);
M.mut.gene_idx = listmap(M.mut.gene,M.gene.name);
idx = find(isnan(M.mut.gene_idx));
if ~isempty(idx)
  fprintf('Omitting %d/%d mutations outside gene set.\n',length(idx),slength(M.mut));
  pr(M.mut,{'gene','patient','classification','type','chr','pos','ref_allele','newbase'},idx(1:min(1e3,length(idx))))
  if length(idx)>1e3, fprintf('...\n'); end
  M.mut = reorder_struct_exclude(M.mut,idx);
end

%memaudit

% CATEGORIES 

% look up context65 if necessary
if isfield(M.mut,'context65')
  M.mut = make_numeric(M.mut,'context65');
else
  fprintf('\nLooking up values in context_and_effect track\n');
  M.mut.context_and_effect = double(M.FWB.context_and_effect.get(M.mut.chr,M.mut.pos));
  M.mut.context_and_effect(M.mut.context_and_effect==-1)=nan;
  M.mut.context65 = nansub(M.context_and_effect.context65,M.mut.context_and_effect);
  % NOTE: M.mut.context_and_effect is not used again anywhere in the code
  %       unless we will eventually add an "enforce_effect" parameter
end

% collapse coverage models, first step: collapse 1885 to 192x5effect
M.cov.CE_terr = collapse_975_to_192x5(M.cov.CE_terr);
M.cov.CE_cov = collapse_975_to_192x5(M.cov.CE_cov);
M.cov.dims_of_gene_effect_cov = 'gene, 192categories, effect (-2 -1 0 1 2)';

if isempty(P.categs_file) || strcmpi(P.categs_file,'none')
  fprintf('Category discovery....\n');
  %M.effect.include_in_category_discovery = grepm('^(syn|mis|non|spl)$',M.effect.name);
  midx = find(nansub(~M.effect.is_indel,M.mut.effect_idx)==1);
  n = hist2d(M.mut.context65(midx),M.mut.newbase_idx(midx),1:65,1:4);
  N = sum(sum(M.cov.CE_cov,3),1)'; %TODO: need to parameterize M.effect.include_in_category_discovery
  N = round(3*N(1:3:end,:));
  N(65)=0;
  Nn = collapse_Nn_65_to_32([N n]);
  PP=[]; PP.max_k = P.number_of_categories_to_discover;
  if ~P.use_new_find_mut_categs
    Ks = find_mut_categs(Nn,PP);
    M.categ = Ks{P.number_of_categories_to_discover};
    M.ncat = slength(M.categ);
  else 
    M.categ = new_find_mut_categs(Nn,PP);
    M.ncat = max(M.categ.categ);
  end
else  % user-supplied categs_file
  fprintf('Using categs_file %s\n',P.categs_file)
  K = load_struct(P.categs_file);
  if all(isfield(K,{'trinuc','newbase','categ'}))
    % NEW STYLE
    K = make_numeric(K,'categ');
    if any(K.categ<1 | isnan(K.categ)), error('bad entries in "categ" field of categs_file'); end
    if slength(K)==96
      K2 = K;
      K2.trinuc = rc(K2.trinuc);
      K2.newbase = rc(K2.newbase);
      K = concat_structs({K,K2});
    end
    if length(unique_combos(K.trinuc,K.newbase))~=192, error('problem with categs_file'); end
    M.categ = K;
    M.ncat = max(K.categ);
  elseif all(isfield(K,{'left','right','from','change','type'}))
    % OLD STYLE
    M.categ = K;
    M.ncat = slength(K);
  else
    error('unknown categs_file format');
  end
end

% collapse coverage models, second step: collapse 192x5effect to ncatx5effect
M.cov.CE_terr = collapse_192_to_categ_set(M.cov.CE_terr,M.categ,2);
M.cov.CE_cov = collapse_192_to_categ_set(M.cov.CE_cov,M.categ,2);
M.cov.dims_of_gene_effect_cov = 'gene, ncat, effect (-2 -1 0 1 2)';

% analyze collapsed categories vs. full-spectrum context_and_effect categories
context975 = map_categories_to_65(M.context_and_effect);
[x, y] = find(sum(assign_65x4_to_categ_set(M.categ), 3) > 0);
h = accumarray(x, y, [], @(x) {x});
h(65) = {NaN};
M.context_and_effect.categ = h(context975);

%also annotate with permutation effects
%NOTE: hardcoded for now.
x = mod(1:slength(M.context_and_effect), 5)'; x(x == 0) = 5;
M.context_and_effect.permutations_effect_idx = M.effect.permutations_effect_idx(x);


%1885-based collapse needs to account for newbase as well.
%TODO: add this back
%M.context_and_effect.Q = collapse_context_and_effect_categories(M.context_and_effect,M.categ);

%memaudit

% assign mutations to categories
fprintf('Assigning mutation categories...\n');
M.mut = rmfield_if_exist(M.mut,{'categ','categ_ignoring_null_categ'});
c = assign_65x4_to_categ_set(M.categ);
c = bsxfun(@times,c,1:M.ncat);
c = squeeze(max(c,[],2));
M.mut.categ_idx = zeros(slength(M.mut),1);
for newbase_idx=1:4
  idx = find(M.mut.newbase_idx==newbase_idx);
  M.mut.categ_idx(idx) = nansub(c(:,newbase_idx),M.mut.context65(idx));
end
idx = find(M.mut.categ_idx==0 & ~M.effect.is_indel(M.mut.effect_idx));
if ~isempty(idx)
  fprintf('Omitting %d/%d mutations of unassigned categ.\n',length(idx),slength(M.mut));
  pr(M.mut,{'gene','patient','classification','type','chr','pos','ref_allele','newbase'},idx(1:min(1e3,length(idx))))
  if length(idx)>1e3, fprintf('...\n'); end
  M.mut = reorder_struct_exclude(M.mut,idx);
end

% geneidx lookup table for speedup
fprintf('Computing geneidx lookup table...\n');
M.mut = sort_struct(M.mut,'gene_idx');
[u ui uj] = unique(M.mut.gene_idx,'first');
h = histc(uj,1:length(u));
M.geneidx = cell(M.ng,1); for i=1:length(u), M.geneidx{u(i)} = as_column(ui(i):ui(i)+h(i)-1);end

%memaudit

% COVARIATES 
fprintf('Loading covariates file.\n');
V = load_struct_specify_string_cols(covariates_file,1);  % gene is string
f = fieldnames(V); if ~strcmp(f{1},'gene'), error('first column of covariates file must be "gene"'); end
M.cvname = f(2:end); M.nv = length(M.cvname);
gidx = listmap(M.gene.name,V.gene);
M.V = nan(M.ng,M.nv);
for vi=1:M.nv, M.V(:,vi) = nansub(V.(M.cvname{vi}),gidx); end

% ____________
% DONE LOADING

fprintf('LOAD finished.  '); toc(loadstart);

%memaudit

% MEMORY REDUCTION   reduces M size from 4Gb to 0.5Gb  (pancan13k.v7.3)
M.mut = rmfield_if_exist(M.mut,{'gene','patient','classification','type','ref_allele','newbase','start','end','context_and_effect'});

%memaudit

%fclose(msgfile);
%
%  function memaudit
%    if P.memory_audit
%      whos
%      disp(datestr(now));
%    end
%  end
%
%  function fprintf(varargin)
%    txt = sprintf(varargin{:});
%    fprintf('%s',txt);
%    fprintf(msgfile,'%s',txt);
%  end

end
