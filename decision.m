function S = decision(mu_post, sig_post, kappa, gamma)
kappa = min(max(kappa,1e-12),1-1e-12);
S = (mu_post + norminv(kappa).*sig_post > gamma);
end