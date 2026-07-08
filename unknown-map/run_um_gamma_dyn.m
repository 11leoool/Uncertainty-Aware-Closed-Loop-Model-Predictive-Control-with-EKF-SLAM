% RUN_UM_GAMMA_DYN  Dynamic-obstacle gamma sweep in the unknown-map framing
% (scenario identical to run_um_dynamic.m, incl. track management).
% gamma = 0 corresponds to cv_fixed.

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
cfg.r_sense = 1.2; cfg.safe_buffer = 0.06;
cfg.o0  = [1.3; 0.15];  cfg.vo0 = [-0.16; 0.16];  cfg.obs_r = 0.15;
cfg.obs_sa2_true = 0.0025;  cfg.obs_sa2_filter = 0.02;   % see run_um_dynamic.m
cfg.track_drop = 3.0;

GAMMAS = [0 0.5 1.0 1.5 2.0 2.5 3.0];

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

nG = numel(GAMMAS);
sweep = struct('gamma',GAMMAS,'coll',zeros(1,nG),'path',zeros(1,nG),'term',zeros(1,nG));
lines = {};
lines{end+1} = sprintf('=== Unknown-map DYNAMIC gamma sweep (M = %d) ===', M);
lines{end+1} = sprintf('%-8s | %-11s | %-9s | %-9s', 'gamma','collision %','path len','term err');
lines{end+1} = repmat('-',1,48);
for gi = 1:nG
    c = cfg; c.gamma = GAMMAS(gi);
    if GAMMAS(gi) == 0, md = 'cv_fixed'; else, md = 'cv_cov'; end
    COL = false(M,1); PL = zeros(M,1); TR = zeros(M,1);
    fprintf('gamma = %.1f ', GAMMAS(gi));
    for k = 1:M
        res = mc_run_trial_um_dyn(md, mpc, c, noises{k});
        COL(k) = res.collided; PL(k) = res.path_len; TR(k) = res.term_true_ref;
    end
    sweep.coll(gi) = 100*mean(COL); sweep.path(gi) = mean(PL); sweep.term(gi) = mean(TR);
    fprintf('-> collisions %4.1f%%, path %.3f, term %.4f\n', ...
        sweep.coll(gi), sweep.path(gi), sweep.term(gi));
    lines{end+1} = sprintf('%6.1f   | %8.1f    | %8.3f | %8.4f', ...
        GAMMAS(gi), sweep.coll(gi), sweep.path(gi), sweep.term(gi)); %#ok<SAGROW>
end
txt = strjoin(lines, newline);
disp(txt);
fid = fopen('um_gamma_dyn_results.txt','w'); fprintf(fid,'%s\n',txt); fclose(fid);
save('um_gamma_dyn_results.mat','sweep','cfg','GAMMAS','M');

f1 = figure('Color','w');
yyaxis left
plot(GAMMAS, sweep.coll, '-o', 'LineWidth',1.6); ylabel('Collision rate (%)');
yyaxis right
plot(GAMMAS, sweep.path, '-s', 'LineWidth',1.6); ylabel('Mean path length (m)');
xlabel('\gamma'); grid on
title('Unknown map, dynamic obstacle: safety-efficiency trade-off');
set(gca,'FontName','Times New Roman','FontSize',12);
saveas(f1,'um_gamma_dyn.png');
fprintf('\nSaved um_gamma_dyn_results.mat/.txt and um_gamma_dyn.png\n');
