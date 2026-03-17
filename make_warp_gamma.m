function W = make_warp_gamma(k,theta)
W.name  = 'Gamma';
W.k = k; W.theta = theta;
W.F     = @(y) gamcdf(y,k,theta);
W.invF  = @(u) gaminv(u,k,theta);
W.phi_inv = @(y) norminv(clamp01(gamcdf(y,k,theta)));
W.meanF = k*theta;
end