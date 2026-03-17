function moments = gaussian_raw_moments(mu,sig,M)
Ng = numel(mu);
s2 = sig.^2;
temp = zeros(Ng,M+1);
temp(:,1)=1;
if M>=1, temp(:,2)=mu; end
if M>=2, temp(:,3)=mu.*temp(:,2)+1*s2.*temp(:,1); end
for m=3:M
    temp(:,m+1)=mu.*temp(:,m)+(m-1)*s2.*temp(:,m-1);
end
moments=temp(:,2:end);
end