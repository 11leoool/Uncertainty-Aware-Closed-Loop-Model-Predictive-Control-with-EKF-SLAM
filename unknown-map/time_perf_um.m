% TIME_PERF_UM  Per-step NMPC solve-time measurement for the unknown-map
% framing: free-space controller (N=10, dt=0.05) on the patrol, and the
% parameterised-obstacle controller (N=14, dt=0.2, slacks) on the sensed
% static obstacle with the covariance-aware margin. 10 trials each.

clear; clc;
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

cfg.lm   = [-0.5 1 -1.5;  0.5 1  0.0];
cfg.wps  = [ 1.5  -1.5   0.0;  1.5   1.0   0.0];
cfg.x0_nom = [0; 0; 0];
cfg.lm_uninf_std = 10;
cfg.var_v = 0.01;  cfg.var_w = 0.01;
cfg.var_d = 0.01;  cfg.var_a = 0.01;
cfg.sim_tim = 45;  cfg.tol = 0.05;  cfg.tol_wp = 0.10;  cfg.wp_timeout = 15;
cfg.r_sense = 1.2; cfg.safe_buffer = 0.06;  cfg.gamma = 2.0;

Mt = 10;
L = size(cfg.lm,2);

% ---- free-space controller ----
mpcF = mc_build_mpc();
maxiterF = round(cfg.sim_tim/mpcF.T);
rng(2024);
tf = [];
for k = 1:Mt
    nz.u = [sqrt(cfg.var_v); sqrt(cfg.var_w)] .* randn(2, maxiterF);
    nz.z = zeros(2, L, maxiterF);
    nz.z(1,:,:) = sqrt(cfg.var_d) * randn(1, L, maxiterF);
    nz.z(2,:,:) = sqrt(cfg.var_a) * randn(1, L, maxiterF);
    r = mc_run_trial_um2('slam', mpcF, cfg, nz);
    tf = [tf, r.solve_t]; %#ok<AGROW>
end
fprintf('Free-space NMPC (N=%d, dt=%.2f): mean %.2f ms, p95 %.2f ms (period %.0f ms)\n', ...
    mpcF.N, mpcF.T, 1000*mean(tf), 1000*prctile(tf,95), 1000*mpcF.T);

% ---- parameterised-obstacle controller ----
optO.N = 14; optO.T = 0.2; optO.rob_diam = 0.3;
optO.v_max = 0.6; optO.omega_max = pi/4; optO.xy_min = -2; optO.xy_max = 2;
mpcO = mc_build_mpc_dyn(optO);
maxiterO = round(cfg.sim_tim/mpcO.T);
cfgO = cfg; cfgO.obs = [0 1.25 0.15];
rng(2024);
to = [];
for k = 1:Mt
    nz.u = [sqrt(cfg.var_v); sqrt(cfg.var_w)] .* randn(2, maxiterO);
    nz.z = zeros(2, L, maxiterO);
    nz.z(1,:,:) = sqrt(cfg.var_d) * randn(1, L, maxiterO);
    nz.z(2,:,:) = sqrt(cfg.var_a) * randn(1, L, maxiterO);
    nz.zo = [sqrt(cfg.var_d); sqrt(cfg.var_a)] .* randn(2, maxiterO);
    r = mc_run_trial_um_obs('slam_cov', mpcO, cfgO, nz);
    to = [to, r.solve_t]; %#ok<AGROW>
end
fprintf('Obstacle NMPC  (N=%d, dt=%.2f): mean %.2f ms, p95 %.2f ms (period %.0f ms)\n', ...
    mpcO.N, mpcO.T, 1000*mean(to), 1000*prctile(to,95), 1000*mpcO.T);

fid = fopen('time_perf_um_results.txt','w');
fprintf(fid, 'Free-space: mean %.2f ms, p95 %.2f ms (period %.0f ms)\n', ...
    1000*mean(tf), 1000*prctile(tf,95), 1000*mpcF.T);
fprintf(fid, 'Obstacle:   mean %.2f ms, p95 %.2f ms (period %.0f ms)\n', ...
    1000*mean(to), 1000*prctile(to,95), 1000*mpcO.T);
fclose(fid);
fprintf('Saved time_perf_um_results.txt\n');
