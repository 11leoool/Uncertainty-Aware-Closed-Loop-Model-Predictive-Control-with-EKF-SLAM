% RUN_UM_DYN_ABLATION  Dynamic matched-mean ablation (review roadmap item 1).
% Phase 1: run cv_cov (gamma = 2) and measure the mean node-0 inflation it
%          actually applies over all trials.
% Phase 2: run cv_fixed with a CONSTANT inflation equal to that mean, on
%          identical seeds. If the size-matched constant is equally safe, the
%          static ablation's conclusion (safety ~ margin size; adaptivity's
%          value = self-tuning) extends to the dynamic case; if not, the
%          adaptivity claim gains direct evidence. Either outcome is reported.
% Scenario identical to run_um_dynamic.m (final configuration).

clear; clc; close all;
addpath('..'); addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

cfg.lm   = [-0.5 1 -1.5;
             0.5 1  0.0];
cfg.wps  = [ 1.5  -1.5   0.0;
             1.5   1.0   0.0];
cfg.x0_nom = [0; 0; 0];
cfg.lm_uninf_std = 10;
cfg.var_v = 0.01;  cfg.var_w = 0.01;
cfg.var_d = 0.01;  cfg.var_a = 0.01;
cfg.sim_tim = 75;  cfg.tol = 0.05;  cfg.tol_wp = 0.10;  cfg.wp_timeout = 30;
cfg.r_sense = 1.2; cfg.safe_buffer = 0.06;  cfg.gamma = 2.0;
cfg.o0  = [1.3; 0.15];  cfg.vo0 = [-0.16; 0.16];  cfg.obs_r = 0.15;
cfg.obs_sa2_true = 0.0025;  cfg.obs_sa2_filter = 0.02;
cfg.track_drop = 3.0;

opt.N = 14; opt.T = 0.2; opt.rob_diam = 0.3;
opt.v_max = 0.6; opt.omega_max = pi/4;
opt.xy_min = -2; opt.xy_max = 2;
fprintf('Building parameterised-obstacle NMPC ...\n');
mpc = mc_build_mpc_dyn(opt);

L = size(cfg.lm,2);
maxiter = round(cfg.sim_tim/mpc.T);
ntot = maxiter + mpc.N + 2;

M = 50;

%% ---- CALIBRATION STREAM (frozen protocol DYN_ABLATION_ESTIMAND.md) ----
% Disjoint seeds rng(81000+t). Estimates c = per-trial mean node-0 inflation
% conditional on an active constraint, equal trial weight, cv_cov arm only.
% The evaluation seeds (rng 2024 family) are NOT touched here.
noises = cell(M,1);
for kk = 1:M
    rng(81000 + kk, 'twister');
    nz.u = [sqrt(cfg.var_v); sqrt(cfg.var_w)] .* randn(2, maxiter);
    nz.z = zeros(2, L, maxiter);
    nz.z(1,:,:) = sqrt(cfg.var_d) * randn(1, L, maxiter);
    nz.z(2,:,:) = sqrt(cfg.var_a) * randn(1, L, maxiter);
    nz.oz   = [sqrt(cfg.var_d); sqrt(cfg.var_a)] .* randn(2, maxiter);
    nz.oacc = sqrt(cfg.obs_sa2_true) * randn(2, ntot);
    noises{kk} = nz;
end

cvals = zeros(M,1); af = zeros(M,1);
for kk = 1:M
    r = mc_run_trial_um_dyn_v2('cv_cov', mpc, cfg, noises{kk});
    cvals(kk) = r.mean_inflation_act; af(kk) = r.active_frac;
end
c = mean(cvals);
fid = fopen('calibration_result.txt','w');
fprintf(fid, 'CALIBRATION (disjoint stream rng(81000+t), t=1..%d)\n', M);
fprintf(fid, 'c = mean over trials of (mean infl | active) = %.4f m\n', c);
fprintf(fid, 'per-trial c: median %.4f  IQR %.4f  min %.4f  max %.4f\n', ...
    median(cvals), iqr(cvals), min(cvals), max(cvals));
fprintf(fid, 'active fraction: mean %.3f  median %.3f\n', mean(af), median(af));
fclose(fid); type calibration_result.txt
save('calibration_result.mat','c','cvals','af','M');
fprintf('\nCALIBRATION_DONE\n');
