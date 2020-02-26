function sigma = gaussian_sigma_from_x_X(x,X)

alpha = x+1;
beta = X-x+1;

sigma = sqrt(alpha.*beta./((alpha+beta).^2.*(alpha+beta+1)));






