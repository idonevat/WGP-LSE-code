function k = kappa_exceed(mu_a,sig_a,mu_b,sig_b,E)
I = double(E>=0);
A = normcdf((E-mu_a)./sig_a);
B = normcdf((E-mu_b)./sig_b);
num = B-I;
den = A+B-2*I;
k = num./den;
k(abs(den)<1e-15)=0.5;
end