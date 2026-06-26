% RUN_GAMMA_SWEEP  Safety-vs-efficiency trade-off of the covariance-aware
% obstacle margin. Sweeps the chance-constraint factor gamma and, for each
% value, runs the cov-aware 'slam' controller over M shared-noise trials,
% recording collision rate, detour (path length) and terminal accuracy.
%
% gamma = 0 reproduces Stage A (fixed margin). odom and oracle are run once as
% gamma-independent reference lines.

clear; clc; close all;
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

% ---------------- scenario (identical to Stage A/B) ----------------
cfg.lm   = [-0.5 1 -1.5; 0.5 1 0.0];
cfg.xs   = [1.5; 1.5; 0];
cfg.x0_nom = [0; 0; 0];
cfg.sigma_init_pos = 0.10; cfg.sigma_init_th = 0.05;
cfg.lm_prior_std   = 0.10;
cfg.var_v = 1e-6; cfg.var_w = 1e-6;
cfg.var_d = 0.01; cfg.var_a = 0.01;
cfg.sim_tim = 30; cfg.tol = 0.05;
cfg.safe_buffer = 0.06;

obstacle = [0.5 0.5 0.15];
opt.N = 14; opt.T = 0.2; opt.rob_diam = 0.3;
opt.v_max = 0.6; opt.omega_max = pi/4;
opt.obs = obstacle; opt.xy_min = -2; opt.xy_max = 2;
fprintf('Building obstacle NMPC ...\n');
mpc = mc_build_mpc_obs(opt);
L = size(cfg.lm,2); maxiter = round(cfg.sim_tim/mpc.T);

gammas = [0 0.5 1.0 1.5 2.0 2.5 3.0];
M = 50;

% ---------------- shared noise (same as Stage A/B) ----------------
rng(2024);
noises = cell(M,1);
for k = 1:M
    nz.offset = [cfg.sigma_init_pos*randn; cfg.sigma_init_pos*randn; cfg.sigma_init_th*randn];
    nz.u = [sqrt(cfg.var_v); sqrt(cfg.var_w)] .* randn(2, maxiter);
    nz.z = zeros(2, L, maxiter);
    nz.z(1,:,:) = sqrt(cfg.var_d) * randn(1, L, maxiter);
    nz.z(2,:,:) = sqrt(cfg.var_a) * randn(1, L, maxiter);
    nz.lm = cfg.lm_prior_std * randn(2, L);
    noises{k} = nz;
end

plen = @(tr) sum(sqrt(sum(diff(tr(1:2,:),1,2).^2,1)));   % path length (m)

% ---------------- reference baselines (gamma-independent) ----------------
c0 = cfg; c0.cov_aware = false;
col_odom = false(M,1); len_oracle = zeros(M,1); len_odom = zeros(M,1);
for k = 1:M
    ro = mc_run_trial_obs('oracle', mpc, c0, noises{k});
    rd = mc_run_trial_obs('odom',   mpc, c0, noises{k});
    len_oracle(k) = plen(ro.traj_true);
    col_odom(k)   = rd.collided;  len_odom(k) = plen(rd.traj_true);
end
odom_rate = 100*mean(col_odom);
oracle_len = mean(len_oracle);

% ---------------- sweep ----------------
G = numel(gammas);
col_rate = zeros(G,1); path_len = zeros(G,1); term_err = zeros(G,1); min_cl = zeros(G,1);
for gi = 1:G
    c = cfg; c.cov_aware = true; c.gamma = gammas(gi);
    COL = false(M,1); PL = zeros(M,1); TR = zeros(M,1); MC = zeros(M,1);
    fprintf('gamma = %.1f ... ', gammas(gi));
    for k = 1:M
        r = mc_run_trial_obs('slam', mpc, c, noises{k});
        COL(k) = r.collided; PL(k) = plen(r.traj_true);
        TR(k) = r.term_true_ref; MC(k) = r.min_clear;
    end
    col_rate(gi) = 100*mean(COL);
    path_len(gi) = mean(PL);
    term_err(gi) = mean(TR);
    min_cl(gi)   = mean(MC);
    fprintf('collision=%4.1f%%  pathlen=%.3f m  term=%.3f m\n', col_rate(gi), path_len(gi), term_err(gi));
end

% ---------------- report ----------------
fprintf('\n=== Gamma sweep (M=%d). Reference: oracle pathlen=%.3f m, odom collision=%.1f%% ===\n', ...
        M, oracle_len, odom_rate);
fprintf('%-6s | %-12s | %-12s | %-12s | %-12s\n','gamma','collision %','pathlen [m]','detour [m]','term err [m]');
fprintf('%s\n', repmat('-',1,64));
for gi = 1:G
    fprintf('%-6.1f | %10.1f   | %10.3f   | %+10.3f   | %10.3f\n', ...
        gammas(gi), col_rate(gi), path_len(gi), path_len(gi)-oracle_len, term_err(gi));
end
fprintf('%s\n', repmat('-',1,64));

save('gamma_sweep.mat','gammas','col_rate','path_len','term_err','min_cl', ...
     'odom_rate','oracle_len','cfg','obstacle','M');
fprintf('Saved gamma_sweep.mat\n');

% ---------------- figure: safety vs efficiency ----------------
f = figure('Color','w','Name','Gamma sweep');
yyaxis left
plot(gammas, col_rate, '-o','linewidth',1.8,'MarkerFaceColor','auto'); hold on
yline(odom_rate,'--','odom','Color',[0.7 0 0],'LineWidth',1.2,'LabelHorizontalAlignment','left');
ylabel('Collision rate (%)'); ylim([-2 max(odom_rate,max(col_rate))+4]);
yyaxis right
plot(gammas, path_len, '-s','linewidth',1.8);
yline(oracle_len,':','oracle path','LineWidth',1.2,'LabelHorizontalAlignment','right');
ylabel('Mean path length (m)');
xlabel('\gamma (chance-constraint factor)');
title('Covariance-aware margin: safety vs. efficiency trade-off');
grid on; set(gca,'FontName','Times New Roman','FontSize',12);
legend({'slam collision rate','odom (reference)','slam path length','oracle path'},'Location','east');
saveas(f,'gamma_sweep.png');
fprintf('Saved gamma_sweep.png\n');
