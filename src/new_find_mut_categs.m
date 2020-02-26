function X = new_find_mut_categs(Nn,P)
% Nn should be a 32x5 table as from collapse_Nn_64_by_strand
%   rows:  32 strand-collapsed categories
%   columns:  [N ->A ->C ->G ->T]
%
% X is list of best category sets discovered, in NEW STYLE


if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'max_k',8);

fprintf('new_find_mut_categs\n');

if size(Nn,1)~=32, error('input must have 32 rows'); end
if size(Nn,2)~=5, error('input must have 5 columns (N A C G T)'); end

Y = generate_categ_context65_names();
Y = reorder_struct(Y,1:32);
Y = rmfield(Y,'num');
Y = parse_in(Y,'name','(.) in (.)_(.)$',{'f','l','r'});
Y.trinuc = regexprep(Y.name,'(.) in (.)_(.)$','$2$1$3');
Y.f = listmap(Y.f,{'A','C','G','T'});
Y.l = listmap(Y.l,{'A','C','G','T'});
Y.r = listmap(Y.r,{'A','C','G','T'});

Y.N = Nn(:,1);

X = cell(4,1);
base = 'ACGT';
for i=1:4
  X{i} = Y;
  X{i}.n = Nn(:,1+i);
  X{i}.t = repmat(i,slength(X{i}),1);
  X{i}.name = regexprep(X{i}.name,'^(.) in (.)_(.)$',['$1->' base(i) ' in $2_$3']);
end
X = concat_structs(X);
X.newbase = nansub({'A','C','G','T'},X.t);
X = order_fields(X,{'name','trinuc','newbase','f','t','l','r','N','n'});
X = reorder_struct_exclude(X,X.f==X.t);

% cluster rates (actually lower bounds on 95% binomial ci)

[X.rate ci] = binofit(X.n,X.N);
X.ci_low = ci(:,1);
X.ci_high = ci(:,2);
X.lograte = log10(X.rate);

%input = X.rate;
input = X.ci_low;
d = pdist(input);
l = linkage(d,'complete');
c = cluster(l,'maxclust',P.max_k);
X.categ = c;

pr(sort_struct(X,'lograte'))



