function outmat = collapse_975_to_192x5(covmat)
%temporary function to map the e15 context+effect track to channels, and collapse out the chromatin states
%should ultimately be replaced with a generic C+E -> C/E converted
%need to respect Mike's order of the 192

%192 ordering (least to most)
% L(B->N)R :: N, R, L, B

%c65 -> ch192
%ch192 -> c65
ch192to65 = NaN(192, 1);
i = 1;
for l = 0:3, for b = 0:3, for n = 0:3, for r = 0:3,
  if b == n, continue; end
  ch192to65(i, 1) = dot([r l b], 4.^(0:2)); i = i + 1; 
end, end, end, end

outmat = NaN(size(covmat, 1), 192, 5);
for i = 1:192,
  outmat(:, i, :) = sum(reshape(covmat(:, (ch192to65(i)*15 + 1):((ch192to65(i) + 1)*15)), size(covmat, 1), [], 3), 3);
end
