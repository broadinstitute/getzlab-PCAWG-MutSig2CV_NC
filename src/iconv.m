function c = iconv(p, nmax)
  %iconv: iterated convolutions
  %inputs: discrete probability dist. p, number of times to iterate nmax
  %outputs: iterated convolution

  if nmax < 1, c = p; return; end

  c = p;
  for i = 1:nmax - 1,
    c = conv(c, p);
  end
end
