%% -------------------------
% Warp constructors
%% -------------------------
function W = make_warp_beta(a,b)
W.name  = 'Beta';
W.a = a; W.b = b;
W.F     = @(y) betacdf(y,a,b);
W.invF  = @(u) betainv(u,a,b);
W.phi_inv = @(y) norminv(clamp01(betacdf(y,a,b)));
W.meanF = a/(a+b);
end
