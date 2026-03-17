%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFECT-MATCH simulations (paper-aligned objectives) vs n
% For each warp (Beta/Gamma/Weibull), we produce ONE figure with 2 panels:
%   (Left)  Moments objective evaluated for the moments-optimal rule
%   (Right) Exceedance objective evaluated for the exceedance-optimal rule
%
% We plot curves for multiple truth/inference lengthscales ell in ellList.
%
% Notes:
% - We compute spatial averages over a uniform grid on [0,1]^2:
%       moemtns based loss: mean_x [ FA(x)*qA(x) + MD(x)*qB(x) ]
%       exceedance based loss:mean_x [ FA(x)*pA(x) + MD(x)*pB(x) ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
rng(11);

R_mc = 100;   % MC repetitions

%% -------------------------
% Observation sizes + MC
%% -------------------------
Nlist=10:10:300;
%% -------------------------
% Grid on [0,1]^2
%% -------------------------
Nx = 50; Ny = 50;
x1 = linspace(0,1,Nx);
x2 = linspace(0,1,Ny);
[X1,X2] = meshgrid(x1,x2);
Xg = [X1(:), X2(:)];
Ng = size(Xg,1);

%% -------------------------
% Lengthscales to test (truth == inference)
%% -------------------------
ellList = [0.05 0.08 0.10 0.15 0.25];   
sig_true = 1.0;
mu0 = 0;

%% -------------------------
% weights and threshold
%% -------------------------
omega = [1 0.5 0.25]; % weights for the moment objective 
E = 0.1;              % threshold for exceedance objective
M = numel(omega);
p = 0.9;  % exceedance quantile for defining lambda
%% -------------------------
% Loss model priors
%% -------------------------
muAlpha = 1; sigAlpha = 1;
muBeta  = 2; sigBeta  = 1;
mu_a_all  = muAlpha*ones(Ng,1);
sig_a_all = sigAlpha*ones(Ng,1);
mu_b_all  = muBeta*ones(Ng,1);
sig_b_all = sigBeta*ones(Ng,1);
%% -------------------------
% Precompute kappas (loss does not depend on warp)
%% -------------------------
k1 = kappa_moments(mu_a_all,sig_a_all,mu_b_all,sig_b_all,omega);   % for WGP-R1
k2 = kappa_exceed (mu_a_all,sig_a_all,mu_b_all,sig_b_all,E);       % for WGP-R2

%% -------------------------
% qA(x)=sum_m ω_m E[α(x)^m],
% qB(x)=sum_m ω_m E[β(x)^m]
%% -------------------------
mA = gaussian_raw_moments(mu_a_all, sig_a_all, M);   % Ng×M
mB = gaussian_raw_moments(mu_b_all, sig_b_all, M);   % Ng×M
qA = mA * omega(:);                                  % Ng×1
qB = mB * omega(:);                                  % Ng×1
%% -------------------------
% pA(x)=P(α(x)>E), 
% pB(x)=P(β(x)>E)
%% -------------------------
pA = 1 - normcdf((E - mu_a_all)./sig_a_all);          % Ng×1
pB = 1 - normcdf((E - mu_b_all)./sig_b_all);          % Ng×1
%% -------------------------
% Warp functions
%% -------------------------
warps = { ...
    make_warp_beta(0.5, 0.5), ...
    make_warp_gamma(0.5, 1.0), ...
    make_warp_weibull(1.0, 5.0) ...
    };

warpNames = { 'Beta ', 'Gamma ', 'Weibull ' };

%% -------------------------
% Fix sensor locations per n (shared across ell + warps for fairness)
%% -------------------------
idxList = cell(numel(Nlist),1);
XnList  = cell(numel(Nlist),1);
for iN = 1:numel(Nlist)
    n = Nlist(iN);
    idxList{iN} = randperm(Ng,n).';
    XnList{iN}  = Xg(idxList{iN},:);
end

%% -------------------------
% Run per warp
%% -------------------------
for iw = 1:numel(warps)
    W = warps{iw};

    % Threshold selection in y-space; corresponding threshold in latent h-space
    lambda    = W.invF(p);            % threshold in y-space
    lambda=0.8; %HARD CODEED 
    gamma_wgp = W.phi_inv(lambda);    % threshold in h-space (for WGP inference)

    % Storage: curves for each ell (rows=|ellList|, cols=|Nlist|)
    J1_R1 = zeros(numel(ellList), numel(Nlist));   % J1 for moments-optimal rule (WGP-R1)
    J2_R2 = zeros(numel(ellList), numel(Nlist));   % J2 for exceedance-optimal rule (WGP-R2)

    % (Optional) MC std for error bars (same dims)
    J1_R1_std = zeros(numel(ellList), numel(Nlist));
    J2_R2_std = zeros(numel(ellList), numel(Nlist));

    for ie = 1:numel(ellList)
        ell_true = ellList(ie);

        % Latent GP truth covariance factor (depends on ell_true)
        Kgg = se_kernel(Xg, Xg, ell_true, sig_true);
        Lgg = chol(Kgg + 1e-10*eye(Ng), 'lower');

        % Inference hyperparams matched to truth
        ell = ell_true;
        sig = sig_true;

        for iN = 1:numel(Nlist)
            n   = Nlist(iN);
            idx = idxList{iN};
            Xn  = XnList{iN};

            j1_mc = zeros(R_mc,1);
            j2_mc = zeros(R_mc,1);

            for r = 1:R_mc
                % ---- Sample latent GP truth ----
                h_true = mu0 + Lgg*randn(Ng,1);

                % ---- Warp to get y_true ----
                u_true = clamp01(normcdf(h_true));
                y_true = W.invF(u_true);

                S_true = (y_true >= lambda);   % Ng×1 logical
                y_obs  = y_true(idx);

                % ---- WGP inference on latent h ----
                z = W.phi_inv(y_obs);
                [mu_post, sig_post] = gp_posterior_noise_free(Xn, z, Xg, mu0, ell, sig);

                % ---- Paper decision rules ----
                S_R1 = decision(mu_post, sig_post, k1, gamma_wgp);  % moments-optimal rule
                S_R2 = decision(mu_post, sig_post, k2, gamma_wgp);  % exceedance-optimal rule

                % ---- PERFECT-MATCH objectives ----
                % FA = (~S_true) & (S_hat),  MD = (S_true) & (~S_hat)
                FA = (~S_true) & (S_R1);  MD = (S_true) & (~S_R1);
                j1_mc(r) = mean( double(FA).*qA + double(MD).*qB );

                FA = (~S_true) & (S_R2);  MD = (S_true) & (~S_R2);
                j2_mc(r) = mean( double(FA).*pA + double(MD).*pB );
            end

            J1_R1(ie,iN) = mean(j1_mc);
            J2_R2(ie,iN) = mean(j2_mc);

            J1_R1_std(ie,iN) = std(j1_mc);
            J2_R2_std(ie,iN) = std(j2_mc);

            fprintf('[%s] ell=%.3f | n=%d | J1(R1)=%.4f | J2(R2)=%.4f\n', ...
                warpNames{iw}, ell_true, n, J1_R1(ie,iN), J2_R2(ie,iN));
        end
    end

    name_WS=['WS_optimal_detection_schemes_' lower(warpNames{iw})];
    save(name_WS)
end

PLOT_OptimalDetectionSchemes
