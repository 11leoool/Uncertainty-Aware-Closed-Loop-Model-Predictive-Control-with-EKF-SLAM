% RUN_UM_GAMMA_ABLATION  Stage 3 of the unknown-map refactor:
%   (a) gamma sweep of the covariance-aware margin (slam_cov), gamma 0..3;
%   (b) matched-mean ablation: a CONSTANT margin equal to slam_cov's mean
%       applied inflation at gamma = 2, on identical seeds (slam feedback).
% Scenario identical to run_um_obstacle.m (patrol, r_sense 1.2, obstacle
% sensed on detection). gamma = 0 reproduces slam_fixed.

clear; clc; close all;
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

% ---------------- scenario (identical to run_um_obstacle.m) ----------------
cfg.lm   = [-0.5 1 -1.5;
             0.5 1  0.0];
cfg.wps  = [ 1.5  -1.5   0.0;
             1.5   1.0   0.0];
cfg.x0_nom = [0; 0; 0];
cfg.obs  = [0 1.25 0.15];

cfg.lm_uninf_std = 10;
cfg.var_v = 0.01;  cfg.var_w = 0.01;
cfg.var_d = 0.01;  cfg.var_a = 0.01;

cfg.sim_tim = 45;
cfg.tol     = 0.05;
cfg.tol_wp  = 0.10;
cfg.wp_timeout = 15;
cfg.r_sense = 1.2;
cfg.safe_buffer = 0.06;

GAMMAS = [0 0.5 1.0 1.5 2.0 2.5 3.0];

% ---------------- controller ----------------
opt.N = 14; opt.T = 0.2; opt.rob_diam = 0.3;
opt.v_max = 0.6; opt.omega_max = pi/4;
opt.xy_min = -2; opt.xy_max = 2;
fprintf('Building parameterised-obstacle NMPC ...\n');
mpc = mc_build_mpc_dyn(opt);

L = size(cfg.lm,2);
maxiter = round(cfg.sim_tim/mpc.T);

% ---------------- shared noise ----------------
M = 50;
rng(2024);
noises = cell(M,1);
for k = 1:M
    nz.u = [sqrt(cfg.var_v); sqrt(cfg.var_w)] .* randn(2, maxiter);
    nz.z = zeros(2, L, maxiter);
    nz.z(1,:,:) = sqrt(cfg.var_d) * randn(1, L, maxiter);
    nz.z(2,:,:) = sqrt(cfg.var_a) * randn(1, L, maxiter);
    nz.zo = [sqrt(cfg.var_d); sqrt(cfg.var_a)] .* randn(2, maxiter);
    noises{k} = nz;
end

% ---------------- (a) gamma sweep ----------------
nG = numel(GAMMAS);
sweep = struct('gamma',GAMMAS,'coll',zeros(1,nG),'path',zeros(1,nG), ...
               'term',zeros(1,nG),'minclear',zeros(1,nG),'meaninfl',zeros(1,nG));
for gi = 1:nG
    c = cfg; c.gamma = GAMMAS(gi);
    if GAMMAS(gi) == 0
        md = 'slam_fixed';                 % gamma=0 == fixed d0
    else
        md = 'slam_cov';
    end
    COL = false(M,1); PL = zeros(M,1); TR = zeros(M,1); MC = zeros(M,1); MI = zeros(M,1);
    fprintf('gamma = %.1f ', GAMMAS(gi));
    for k = 1:M
        res = mc_run_trial_um_obs(md, mpc, c, noises{k});
        COL(k) = res.collided; PL(k) = res.path_len;
        TR(k) = res.term_true_ref; MC(k) = res.min_clear;
        MI(k) = res.mean_inflation;
    end
    sweep.coll(gi) = 100*mean(COL);  sweep.path(gi) = mean(PL);
    sweep.term(gi) = mean(TR);       sweep.minclear(gi) = min(MC);
    sweep.meaninfl(gi) = mean(MI);
    fprintf('-> collisions %4.1f%%, path %.3f, term %.4f\n', ...
        sweep.coll(gi), sweep.path(gi), sweep.term(gi));
end

% ---------------- (b) matched-mean ablation ----------------
gi2 = find(GAMMAS==2.0);
marg = sweep.meaninfl(gi2);          % mean applied inflation at gamma = 2
fprintf('\nMatched fixed margin (slam_cov gamma=2 mean inflation): %.4f m\n', marg);
c = cfg; c.gamma = 0; c.fixed_extra = marg;
COL = false(M,1); PL = zeros(M,1); TR = zeros(M,1); MC = zeros(M,1);
fprintf('fixed-matched ');
for k = 1:M
    res = mc_run_trial_um_obs('slam_fixed', mpc, c, noises{k});
    COL(k) = res.collided; PL(k) = res.path_len;
    TR(k) = res.term_true_ref; MC(k) = res.min_clear;
end
abl.coll = 100*mean(COL); abl.path = mean(PL);
abl.term = mean(TR); abl.minclear = min(MC); abl.marg = marg;
fprintf('-> collisions %4.1f%%, path %.3f\n', abl.coll, abl.path);

% ---------------- report ----------------
lines = {};
lines{end+1} = sprintf('=== Unknown-map gamma sweep + ablation (M = %d) ===', M);
lines{end+1} = sprintf('%-8s | %-11s | %-9s | %-9s | %-10s | %-10s', ...
    'gamma','collision %','path len','term err','min clear','mean infl');
lines{end+1} = repmat('-',1,72);
for gi = 1:nG
    lines{end+1} = sprintf('%6.1f   | %8.1f    | %8.3f | %8.4f | %+8.4f  | %8.4f', ...
        GAMMAS(gi), sweep.coll(gi), sweep.path(gi), sweep.term(gi), ...
        sweep.minclear(gi), sweep.meaninfl(gi)); %#ok<SAGROW>
end
lines{end+1} = '';
lines{end+1} = sprintf('fixed-matched (%.3f m const): collisions %.1f%%, path %.3f, term %.4f, min clear %+.4f', ...
    abl.marg, abl.coll, abl.path, abl.term, abl.minclear);
lines{end+1} = sprintf('cov gamma=2 (reference)     : collisions %.1f%%, path %.3f, term %.4f, min clear %+.4f', ...
    sweep.coll(gi2), sweep.path(gi2), sweep.term(gi2), sweep.minclear(gi2));
txt = strjoin(lines, newline);
disp(txt);
fid = fopen('um_gamma_ablation_results.txt','w'); fprintf(fid,'%s\n',txt); fclose(fid);
save('um_gamma_ablation_results.mat','sweep','abl','cfg','GAMMAS','M');

% ---------------- figure ----------------
f1 = figure('Color','w');
yyaxis left
plot(GAMMAS, sweep.coll, '-o', 'LineWidth',1.6); ylabel('Collision rate (%)');
yyaxis right
plot(GAMMAS, sweep.path, '-s', 'LineWidth',1.6); ylabel('Mean path length (m)');
xlabel('\gamma'); grid on
title('Unknown map: safety-efficiency trade-off of the covariance-aware margin');
set(gca,'FontName','Times New Roman','FontSize',12);
saveas(f1,'um_gamma_sweep.png');
fprintf('\nSaved um_gamma_ablation_results.mat/.txt and um_gamma_sweep.png\n');
