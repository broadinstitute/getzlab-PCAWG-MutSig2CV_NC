function pCF = calc_pCF(M,ng,np,P)
% calc_af_score(<mutation struct M>, <ng>, <np>, <parameters P>)
%
% returns [<pCF>]
%
% 2014-11-24//Julian Hess 

if nargin < 1, error('Not enough arguments.'); end

if ~exist('P'), P = []; end;
P = impose_default_value(P, 'min_effect_size', 1.00);
P = impose_default_value(P, 'transform_pmf', 0);
P = impose_default_value(P, 'pmf_trans_fun', 'uniform');
P = impose_default_value(P, 'pmf_trans_fun_p', 50);
P = impose_default_value(P, 'blur_scores', 0);
P = impose_default_value(P, 'blur_std', 0.16);

demand_fields(M,{'gene_idx','pat_idx','i_tumor_f'});
M = make_numeric(M,{'gene_idx','pat_idx','i_tumor_f'});
M = sort_struct(M, {'pat_idx', 'gene_idx'});

disp(sprintf('Creating score matrix (%d patients) ...', np));
Gm = NaN(ng,np);
for pi=1:np
  midx = find(M.pat_idx==pi);
  [ug ugi ugj] = unique(M.gene_idx(midx));
  for gi=1:length(ugi)
    gmidx = midx(ugj==gi);
    if ~isempty(gmidx)
      Gm(ug(gi), pi) = 100*nanmean(M.i_tumor_f(gmidx));
    end
  end
  if ~mod(pi, 100), fprintf('%d/%d  ', pi,np); end
end, fprintf('\n');


if P.blur_scores,
  rng(now, 'twister');
  pidx = Gm > 0;

  for i = 1:size(Gm, 2)
    idx = find(Gm(:, i) > 0);
    Gm(idx, i) = Gm(idx, i) + randi(P.blur_std*[-100, 100], length(idx), 1)/100;
  end

  Gm(pidx) = Gm(pidx) + abs(min(min(Gm)));
end

fprintf('Calculating patient PMFs ...\n');
Gm = round(Gm);
Gmt_tmp = round(nansum(Gm, 2));

%make PMFs
m = max(max(Gm));
Gh = NaN(size(Gm, 2), m + 1); %101 -> m + 1
for i = 1:size(Gm, 2),
  Gh(i, :) = [0; histc(Gm(:, i), 1:m)]; %100 -> m
end

%blur distribution (if requested)
if P.transform_pmf,
  disp('Transforming patient PMFs ...');

  p = P.pmf_trans_fun_p;
  if strcmp(P.pmf_trans_fun, 'gaussian'),
    ker = (1/sqrt(2*pi*p^2))*exp(-(-15:0.5:15).^2/(2*p^2));
  elseif strcmp(P.pmf_trans_fun, 'uniform'),
    ker = ones(1, bitor(p, 1));
  end

  kerw = (length(ker) - 1)/2;
  Gh = padarray(Gh, [0 kerw], 0);
  for i = 1:size(Gh, 1),
    Gh(i, :) = pointwise_convolve(Gh(i, (kerw + 1):(end - kerw)), ker);
    Gh(i, :) = Gh(i, :)/sum(Gh(i, :)); %TODO: normalize kernels to avoid this step?
  end

  %compute distribution of sums 
  disp('Calculating gene score distributions ...');
  midx = find(nansum(Gm, 2));
  Sd = cell(length(midx), 1); %sparse matrix buffer

  for ii = 1:length(midx)
    i = midx(ii);

    cvec = 1; 
    for j = find(Gm(i, :) > 0),
      v = histc(Gm(i, j), 0:100);

      p = pointwise_convolve(v, ker); p = p/sum(p); %TODO: normalized kernel to avoid this problem?
      cvec = conv(p, cvec);
    end
    cvec = cvec/sum(cvec);
    idx = cvec > 0;

    Sd{ii} = [i*ones(nnz(idx), 1) find(idx)' cvec(idx)']; 
  end 
  Gmt = cat(1, Sd{:});
else
  %normalize untransformed PMFs
  for i = 1:size(Gh, 1),
    Gh(i, :) = Gh(i, :)/sum(Gh(i, :));
  end

  %distribution of sums (kronecker)
  Gmt = [find(Gmt_tmp > 0) Gmt_tmp(Gmt_tmp > 0) ones(nnz(Gmt_tmp), 1)];
end

%generate sparse matrix index
[Gmt_idx(:, 1), Gmt_idx(:, 2)] = unique(Gmt(:, 1));
Gmt_idx = [Gmt_idx; [Gmt_idx(end, 1) size(Gmt, 1)]];

%do convolutions
disp('Calculating gene p-values ...');
pdist = cell(size(Gmt_idx, 1) - 1, 1);
pconvm = NaN(size(Gmt_idx, 1) - 1, 1);

iimax = size(Gmt_idx, 1) - 1;
for ii = 1:iimax,
  if ~mod(ii,1000), fprintf('%d/%d ',ii,iimax); end

  i = Gmt_idx(ii, 1);

  cvec = 1;
  for j = find(Gm(i, :) > 0),
    cvec = conv(Gh(j, :), cvec);
  end

  idx = Gmt_idx(ii, 2):(Gmt_idx(ii + 1, 2) - 1);

  sums = [Gmt(idx(1:end), 2); length(cvec)];

  tmp = NaN(length(idx), 1);
  surv = sum(cvec(ceil((sums(1))/P.min_effect_size):end));
  cdf = 0;

  for v = 1:(length(sums) - 1),
    tmp(v) = surv - cdf;
    cdf = cdf + sum(cvec(ceil((sums(v))/P.min_effect_size):(sums(v + 1) - 1)));
  end

  pdist{ii} = [tmp Gmt(idx, 3)]; %save the distribution of p-values if we want it later 
  pconvm(ii) = dot(tmp, Gmt(idx, 3)); %mean pvalue

end, fprintf('\n');

%save output
pCF = ones(ng, 1);
pCF(Gmt_idx(1:(end - 1), 1)) = pconvm;




