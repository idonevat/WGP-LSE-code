function k = kappa_moments(mu_a,sig_a,mu_b,sig_b,omega)
M = numel(omega);
ma = gaussian_raw_moments(mu_a,sig_a,M);
mb = gaussian_raw_moments(mu_b,sig_b,M);
qa = ma*omega(:);
qb = mb*omega(:);
den = qa+qb;
k = qb./den;
k(abs(den)<1e-15)=0.5;
end