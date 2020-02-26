function v = pointwise_convolve(vec, ker)
  if ~bitand(ker, 1), error('Kernel length must be odd.'); end

  w = (length(ker) - 1)/2;

  v = zeros(1, length(vec) + 2*w);
  for j = find(vec),
    idx = ((j - w):(j + w)) + w;
    v(idx) = v(idx) + vec(j)*ker;
  end
end 
