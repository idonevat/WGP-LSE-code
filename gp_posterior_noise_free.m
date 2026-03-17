function [mu_post, sig_post] = gp_posterior_noise_free(Xn, yn, Xg, mu0, ell, sig)
% Exact GP posterior mean/variance at Xg given noise-free observations (Xn, yn)
n = size(Xn,1);
Knn = se_kernel(Xn,Xn,ell,sig) + 1e-10*eye(n);
Lnn = chol(Knn,'lower');
alpha = Lnn'\(Lnn\(yn - mu0));
Kng = se_kernel(Xn,Xg,ell,sig);

mu_post = mu0 + Kng'*alpha;

V = Lnn\Kng;
sig2 = sig^2 - sum(V.^2,1)';
sig2 = max(sig2,0);
sig_post = sqrt(sig2);
end