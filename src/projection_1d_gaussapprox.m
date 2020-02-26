function [p mu sigma] = projection_1d_gaussapprox(Sdeg,Pdeg,score_obs,logspace_fit)

if ~exist('logspace_fit','var'), logspace_fit = false; end

if ~logspace_fit

  mu_pat = sum(Sdeg.*Pdeg,2); var_pat = sum(Pdeg.*(bsxfun(@minus,Sdeg,mu_pat).^2),2);
  mu = sum(mu_pat); var = sum(var_pat); sigma = sqrt(var);
  p = 1-normcdf(score_obs,mu,sigma);

else

  Sdeg = max(Sdeg,0.1);
  mu_pat = sum(Sdeg.*Pdeg,2);
  mu = sum(mu_pat);
  var_pat = sum(Pdeg.*(bsxfun(@minus,log(Sdeg),log(mu_pat)).^2),2);
  var = sum(var_pat); sigma = sqrt(var);
  p = 1-normcdf(log(score_obs),log(mu),sigma);

end
