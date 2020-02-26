function [p am] = pCLmn(x0, nmax, L)
  %pCLmn: calculates pCL using fast multinomial method.
  %inputs: plurality threshold x0, total number of mutations nmax, gene length L
  %output: p-value

  if x0 > nmax, p = NaN; am = NaN; return; end

  a = zeros(nmax + 1, 1); %memoization table for total sum, +1 indexed.  zeros because log(1) = 0
  for n = 1:nmax,
    u = min([n x0 - 1]); %do the sum up to x0 - 1 max muts
    s = zeros(u, 1);
    for i = 1:u,
      s(i) = -sum(log(1:i)) + a(n - i + 1) + log(i*(L + 1) - n);
    end
    m = max(s);
    if isempty(m), s = 0; m = 0; end

    a(n + 1) = -log(n) + m + log(sum(exp(s - m))); %recall that a(i) is actually log(a(i))
  end 

  am = a(nmax + 1);
  p = 1 - exp(sum(log(1:nmax)) - nmax*log(L) + a(nmax + 1));
end
