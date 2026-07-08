% RUN_UM_M100_DYN  M = 100 robustness re-run of the DYNAMIC study only, with
% the corrected mission-executive parameters (wp_timeout = 30 s, sim_tim = 75 s;
% see run_um_dynamic.m header) and the waypoint-achievement metric. Replaces
% the dynamic rows of the earlier um_m100 results; the free-space/static/
% ablation parts of um_m100_results.txt are unaffected by this change.

clear; clc; close all;
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

M = 100;

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
cfg.obs_sa2_true = 0.0025;  cfg.obs_sa2_filter = 0.02;   % see run_um_dynamic.m
cfg.track_drop = 3.0;

opt.N = 14; opt.T = 0.2; opt.rob_diam = 0.3;
opt.v_max = 0.6; opt.omega_max = pi/4;
opt.xy_min = -2; opt.xy_max = 2;
fprintf('Building parameterised-obstacle NMPC ...\n');
mpc = mc_build_mpc_dyn(opt);

L = size(cfg.lm,2);
maxiter = round(cfg.sim_tim/mpc.T);
ntot = maxiter + mpc.N + 2;

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

modes = {'oracle','static','cv_fixed','cv_cov'};
lines = {};
lines{end+1} = sprintf('=== Dynamic M = %d (wp_timeout 30 s, sim_tim 75 s) ===', M);
for mi = 1:numel(modes)
    mode = modes{mi};
    COL = false(M,1); PL = zeros(M,1); TR = zeros(M,1); WV = zeros(M,1);
    fprintf('Running mode: %-9s ', mode);
    for k = 1:M
        r = mc_run_trial_um_dyn(mode, mpc, cfg, noises{k});
        COL(k) = r.collided; PL(k) = r.path_len; TR(k) = r.term_true_ref;
        wv = 0;
        for w = 1:size(cfg.wps,2)
            if min(vecnorm(r.traj(1:2,:) - cfg.wps(:,w))) < 0.15, wv = wv+1; end
        end
        WV(k) = wv;
    end
    [lo,hi] = wilson(sum(COL), M);
    lines{end+1} = sprintf('%-9s coll %5.1f%% [%4.1f, %4.1f] | path %.3f | term %.4f | wps %.2f/3', ...
        mode, 100*mean(COL), 100*lo, 100*hi, mean(PL), mean(TR), mean(WV)); %#ok<SAGROW>
    fprintf('done (collisions: %d/%d, wps %.2f/3).\n', sum(COL), M, mean(WV));
end
txt = strjoin(lines, newline);
disp(txt);
fid = fopen('um_m100_dyn_results.txt','w'); fprintf(fid,'%s\n',txt); fclose(fid);
save('um_m100_dyn_results.mat','lines','cfg','M');
fprintf('Saved um_m100_dyn_results.txt / .mat\n');

function [lo,hi] = wilson(k, n)
z = 1.96; p = k/n;
den = 1 + z^2/n;
ctr = (p + z^2/(2*n)) / den;
hw  = z*sqrt(p*(1-p)/n + z^2/(4*n^2)) / den;
lo = max(0, ctr-hw); hi = min(1, ctr+hw);
end
