% RUN_UM_M100  Consolidated M = 100 robustness re-run of the unknown-map
% studies (free-space patrol, static sensed obstacle, matched-fixed ablation,
% dynamic obstacle). Same rng(2024) stream conventions: the first 50 noise
% structs coincide with the M = 50 primaries. Reports collision rates with
% Wilson 95% CIs and terminal errors.

clear; clc; close all;
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

M = 100;
OUT = {};

% ================= shared scenario pieces =================
base.lm   = [-0.5 1 -1.5;
              0.5 1  0.0];
base.wps  = [ 1.5  -1.5   0.0;
              1.5   1.0   0.0];
base.x0_nom = [0; 0; 0];
base.lm_uninf_std = 10;
base.var_v = 0.01;  base.var_w = 0.01;
base.var_d = 0.01;  base.var_a = 0.01;
base.sim_tim = 45;  base.tol = 0.05;  base.tol_wp = 0.10;  base.wp_timeout = 15;
base.r_sense = 1.2; base.safe_buffer = 0.06;  base.gamma = 2.0;
L = size(base.lm,2);

% ================= (1) free-space patrol =================
fprintf('=== (1) free-space patrol, M = %d ===\n', M);
mpcF = mc_build_mpc();
maxiterF = round(base.sim_tim/mpcF.T);
rng(2024);
noisesF = cell(M,1);
for k = 1:M
    nz.u = [sqrt(base.var_v); sqrt(base.var_w)] .* randn(2, maxiterF);
    nz.z = zeros(2, L, maxiterF);
    nz.z(1,:,:) = sqrt(base.var_d) * randn(1, L, maxiterF);
    nz.z(2,:,:) = sqrt(base.var_a) * randn(1, L, maxiterF);
    noisesF{k} = nz;
end
OUT{end+1} = sprintf('--- Free-space patrol (M = %d) ---', M);
for md = {'oracle','odom','slam'}
    TR = zeros(M,1); LE = zeros(M,1);
    for k = 1:M
        r = mc_run_trial_um2(md{1}, mpcF, base, noisesF{k});
        TR(k) = r.term_true_ref; LE(k) = r.mean_true_est;
    end
    OUT{end+1} = sprintf('%-8s term %.4f +/- %.4f | locerr %.4f', ...
        md{1}, mean(TR), 1.96*std(TR)/sqrt(M), mean(LE));
    fprintf('%s\n', OUT{end});
end

% ================= (2) static sensed obstacle =================
fprintf('=== (2) static obstacle, M = %d ===\n', M);
optO.N = 14; optO.T = 0.2; optO.rob_diam = 0.3;
optO.v_max = 0.6; optO.omega_max = pi/4; optO.xy_min = -2; optO.xy_max = 2;
mpcO = mc_build_mpc_dyn(optO);
maxiterO = round(base.sim_tim/mpcO.T);
cfgO = base; cfgO.obs = [0 1.25 0.15];
rng(2024);
noisesO = cell(M,1);
for k = 1:M
    nz.u = [sqrt(base.var_v); sqrt(base.var_w)] .* randn(2, maxiterO);
    nz.z = zeros(2, L, maxiterO);
    nz.z(1,:,:) = sqrt(base.var_d) * randn(1, L, maxiterO);
    nz.z(2,:,:) = sqrt(base.var_a) * randn(1, L, maxiterO);
    nz.zo = [sqrt(base.var_d); sqrt(base.var_a)] .* randn(2, maxiterO);
    noisesO{k} = nz;
end
OUT{end+1} = sprintf('--- Static sensed obstacle (M = %d) ---', M);
MEANINFL = 0;
for md = {'oracle','odom_fixed','slam_fixed','slam_cov'}
    COL = false(M,1); MI = zeros(M,1); PL = zeros(M,1);
    for k = 1:M
        r = mc_run_trial_um_obs(md{1}, mpcO, cfgO, noisesO{k});
        COL(k) = r.collided; MI(k) = r.mean_inflation; PL(k) = r.path_len;
    end
    [lo,hi] = wil(sum(COL), M);
    OUT{end+1} = sprintf('%-11s coll %5.1f%% [%4.1f, %4.1f] | path %.3f', ...
        md{1}, 100*mean(COL), 100*lo, 100*hi, mean(PL));
    fprintf('%s\n', OUT{end});
    if strcmp(md{1},'slam_cov'), MEANINFL = mean(MI); end
end

% ================= (3) matched-fixed ablation =================
fprintf('=== (3) matched-fixed ablation, M = %d ===\n', M);
cfgA = cfgO; cfgA.fixed_extra = MEANINFL;
COL = false(M,1); PL = zeros(M,1);
for k = 1:M
    r = mc_run_trial_um_obs('slam_fixed', mpcO, cfgA, noisesO{k});
    COL(k) = r.collided; PL(k) = r.path_len;
end
[lo,hi] = wil(sum(COL), M);
OUT{end+1} = sprintf('--- Ablation: fixed-matched %.3f m: coll %.1f%% [%4.1f, %4.1f] | path %.3f', ...
    MEANINFL, 100*mean(COL), 100*lo, 100*hi, mean(PL));
fprintf('%s\n', OUT{end});

% ================= (4) dynamic obstacle =================
fprintf('=== (4) dynamic obstacle, M = %d ===\n', M);
cfgD = base;
cfgD.o0 = [1.3; 0.15]; cfgD.vo0 = [-0.16; 0.16]; cfgD.obs_r = 0.15;
cfgD.obs_sa2_true = 0.0225; cfgD.obs_sa2_filter = 0.02; cfgD.track_drop = 3.0;
ntot = maxiterO + mpcO.N + 2;
rng(2024);
noisesD = cell(M,1);
for k = 1:M
    nz.u = [sqrt(base.var_v); sqrt(base.var_w)] .* randn(2, maxiterO);
    nz.z = zeros(2, L, maxiterO);
    nz.z(1,:,:) = sqrt(base.var_d) * randn(1, L, maxiterO);
    nz.z(2,:,:) = sqrt(base.var_a) * randn(1, L, maxiterO);
    nz.oz   = [sqrt(base.var_d); sqrt(base.var_a)] .* randn(2, maxiterO);
    nz.oacc = sqrt(cfgD.obs_sa2_true) * randn(2, ntot);
    noisesD{k} = nz;
end
OUT{end+1} = sprintf('--- Dynamic obstacle (M = %d) ---', M);
for md = {'oracle','static','cv_fixed','cv_cov'}
    COL = false(M,1); PL = zeros(M,1); TR = zeros(M,1);
    for k = 1:M
        r = mc_run_trial_um_dyn(md{1}, mpcO, cfgD, noisesD{k});
        COL(k) = r.collided; PL(k) = r.path_len; TR(k) = r.term_true_ref;
    end
    [lo,hi] = wil(sum(COL), M);
    OUT{end+1} = sprintf('%-9s coll %5.1f%% [%4.1f, %4.1f] | path %.3f | term %.4f', ...
        md{1}, 100*mean(COL), 100*lo, 100*hi, mean(PL), mean(TR));
    fprintf('%s\n', OUT{end});
end

% ================= save =================
txt = strjoin(OUT, newline);
fid = fopen('um_m100_results.txt','w'); fprintf(fid,'%s\n',txt); fclose(fid);
save('um_m100_results.mat','OUT','M');
fprintf('\nSaved um_m100_results.txt / .mat\n');

% =========================== local functions ===============================
function [lo,hi] = wil(k, n)
z = 1.96; p = k/n;
den = 1 + z^2/n;
ctr = (p + z^2/(2*n)) / den;
hw  = z*sqrt(p*(1-p)/n + z^2/(4*n^2)) / den;
lo = max(0, ctr-hw); hi = min(1, ctr+hw);
end
