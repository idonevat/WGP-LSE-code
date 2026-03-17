function W = make_warp_weibull(A, kshape)
% MATLAB's wblcdf/wblinv use scale A and shape B=kshape
W.name  = 'Weibull';
W.A = A; W.k = kshape;
W.F     = @(y) wblcdf(y,A,kshape);
W.invF  = @(u) wblinv(u,A,kshape);
W.phi_inv = @(y) norminv(clamp01(wblcdf(y,A,kshape)));
W.meanF = A * gamma(1 + 1/kshape);
end