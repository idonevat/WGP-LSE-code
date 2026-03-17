function K = se_kernel(X1,X2,ell,sigma)
D2 = pdist2(X1,X2,'squaredeuclidean');
K = (sigma^2)*exp(-0.5*D2/ell^2);
end
























