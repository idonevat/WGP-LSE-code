function u = clamp01(u)
u = min(max(u,1e-12), 1-1e-12);
end