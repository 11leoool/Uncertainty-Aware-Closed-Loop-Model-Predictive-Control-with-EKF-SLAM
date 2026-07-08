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
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

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
rng(2024);
noises = cell(M,1);
for k = 1:M
    nz.u = [sqrt(cfg.var_v); sqrt(cfg.var_w)] .* randn(2, maxiter);
    nz.z = zeros(2, L, maxiter);
    nz.z(1,:,:) = sqrt(cfg.var_d) * randn(1, L, maxiter);
    nz.z(2,:,:) = sqrt(cfg.var_a) * randn(1, L, maxiter);
    nz.oz   = [sqrt(cfg.var_d); sqrt(cfg.var_a)] .* randn(2, maxiter);
    nz.oacc = sqrt(cfg.obs_sa2_true) * randn(2, ntot);
    noises{k} = nz;
end

lines = {};
lines{end+1} = sprintf('=== Dynamic matched-mean ablation (M = %d) ===', M);

% ---- phase 1: cv_cov, measure mean applied inflation ----
[cov_stats, MI] = run_cell('cv_cov', 0, mpc, cfg, noises, M);
marg = mean(MI);
lines{end+1} = sprintf('cv_cov (gamma=2)      : %s | mean applied inflation %.4f m', cov_stats, marg);
fprintf('%s\n', lines{end});

% ---- phase 2: fixed-matched constant, same seeds ----
[fix_stats, ~] = run_cell('cv_fixed', marg, mpc, cfg, noises, M);
lines{end+1} = sprintf('fixed-matched (%.3f m): %s', marg, fix_stats);
fprintf('%s\n', lines{end});

txt = strjoin(lines, newline);
fid = fopen('um_dyn_ablation_results.txt','w'); fprintf(fid,'%s\n',txt); fclose(fid);
save('um_dyn_ablation_results.mat','lines','marg','cfg','M');
fprintf('Saved um_dyn_ablation_results.txt / .mat\n');

% =========================== local functions ===============================
function [s, MI] = run_cell(md, fixed_extra, mpc, cfg, noises, M)
c = cfg;
if fixed_extra > 0, c.fixed_extra = fixed_extra; end
COL = false(M,1); PL = zeros(M,1); TR = zeros(M,1); WV = zeros(M,1);
MC = zeros(M,1); MI = zeros(M,1);
fprintf('Running %-9s ', md);
for k = 1:M
    r = mc_run_trial_um_dyn(md, mpc, c, noises{k});
    COL(k) = r.collided; PL(k) = r.path_len; TR(k) = r.term_true_ref;
    MC(k) = r.min_clear; MI(k) = r.mean_inflation;
    wv = 0;
    for w = 1:size(cfg.wps,2)
        if min(vecnorm(r.traj(1:2,:) - cfg.wps(:,w))) < 0.15, wv = wv+1; end
    end
    WV(k) = wv;
end
[lo,hi] = wilson(sum(COL), M);
s = sprintf('coll %4.1f%% [%4.1f, %4.1f] | min clear %+7.4f | path %6.3f | term %.4f | wps %.2f/3', ...
    100*mean(COL), 100*lo, 100*hi, min(MC), mean(PL), mean(TR), mean(WV));
fprintf('done.\n');
end

function [lo,hi] = wilson(k, n)
z = 1.96; p = k/n;
den = 1 + z^2/n;
ctr = (p + z^2/(2*n)) / den;
hw  = z*sqrt(p*(1-p)/n + z^2/(4*n^2)) / den;
lo = max(0, ctr-hw); hi = min(1, ctr+hw);
end
