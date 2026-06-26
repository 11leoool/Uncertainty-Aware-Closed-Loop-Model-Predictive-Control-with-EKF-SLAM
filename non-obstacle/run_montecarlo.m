% RUN_MONTECARLO  Monte-Carlo comparison of MPC feedback sources:
%   oracle (true pose)  vs  odom (dead reckoning)  vs  slam (EKF-SLAM).
%
% Each trial uses a random true start offset and a fixed noise realisation that
% is SHARED across the three modes, so differences are due to the estimator
% only. Reports mean +/- 95% CI of the terminal true-vs-reference error
% (how close the REAL robot gets to the goal) and the localisation error.

clear; clc; close all;
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

% ---------------- configuration ----------------
cfg.lm   = [-0.5 1 -1.5;        % landmark x
             0.5 1  0.0];       % landmark y
cfg.xs   = [1.5; 1.5; 0];       % reference posture
cfg.x0_nom = [0; 0; 0];         % nominal start (estimator initialised here)

cfg.sigma_init_pos = 0.10;      % m   true start offset / initial pose uncertainty
cfg.sigma_init_th  = 0.05;      % rad
cfg.lm_prior_std   = 0.10;      % m   prior std of the surveyed landmark map

cfg.var_v = 1e-6;  cfg.var_w = 1e-6;   % control/actuation noise variances
cfg.var_d = 0.01;  cfg.var_a = 0.01;   % measurement noise variances (range,bearing)

cfg.sim_tim = 20;               % s, max sim time
cfg.tol     = 0.05;             % m, stop when feedback within tol of goal

M     = 50;                     % Monte-Carlo trials
modes = {'oracle','odom','slam'};

% ---------------- build controller once ----------------
fprintf('Building NMPC ...\n');
mpc = mc_build_mpc();
L = size(cfg.lm,2);
maxiter = round(cfg.sim_tim/mpc.T);

% ---------------- pre-draw noise (shared across modes) ----------------
rng(2024);                      % master seed -> fully reproducible
noises = cell(M,1);
for k = 1:M
    nz.offset = [cfg.sigma_init_pos*randn; cfg.sigma_init_pos*randn; cfg.sigma_init_th*randn];
    nz.u = [sqrt(cfg.var_v); sqrt(cfg.var_w)] .* randn(2, maxiter);
    nz.z = zeros(2, L, maxiter);
    nz.z(1,:,:) = sqrt(cfg.var_d) * randn(1, L, maxiter);
    nz.z(2,:,:) = sqrt(cfg.var_a) * randn(1, L, maxiter);
    nz.lm = cfg.lm_prior_std * randn(2, L);     % surveyed-map prior error
    noises{k} = nz;
end

% ---------------- run ----------------
results = struct();
for mi = 1:numel(modes)
    mode = modes{mi};
    TR = zeros(M,1); RR = zeros(M,1); LE = zeros(M,1); ST = zeros(M,1);
    sample_traj = [];
    fprintf('Running mode: %-7s ', mode);
    for k = 1:M
        res = mc_run_trial(mode, mpc, cfg, noises{k});
        TR(k) = res.term_true_ref;
        RR(k) = res.rmse_true_ref;
        LE(k) = res.mean_true_est;
        ST(k) = res.steps;
        if k == 1, sample_traj = res; end
    end
    results.(mode).term_true_ref = TR;
    results.(mode).rmse_true_ref = RR;
    results.(mode).loc_err       = LE;
    results.(mode).steps         = ST;
    results.(mode).sample        = sample_traj;
    fprintf('done.\n');
end

% ---------------- report ----------------
ci = @(v) 1.96*std(v)/sqrt(numel(v));
fprintf('\n=== Monte-Carlo results (M = %d) ===\n', M);
fprintf('%-8s | %-22s | %-22s | %-8s\n', 'mode', ...
        'term true-vs-ref [m]', 'mean loc err t-e [m]', 'steps');
fprintf('%s\n', repmat('-',1,72));
for mi = 1:numel(modes)
    m  = modes{mi};
    TR = results.(m).term_true_ref;
    LE = results.(m).loc_err;
    fprintf('%-8s | %6.4f +/- %6.4f      | %6.4f +/- %6.4f      | %5.1f\n', ...
            m, mean(TR), ci(TR), mean(LE), ci(LE), mean(results.(m).steps));
end
fprintf('%s\n', repmat('-',1,72));

save('mc_results.mat','results','cfg','M');
fprintf('Saved mc_results.mat\n');

% ---------------- figures ----------------
% (1) terminal true-vs-reference error: bar + 95%% CI
f1 = figure('Color','w','Name','Terminal true-vs-ref error');
means = cellfun(@(m) mean(results.(m).term_true_ref), modes);
cis   = cellfun(@(m) ci(results.(m).term_true_ref),   modes);
b = bar(means,0.6); hold on
errorbar(1:numel(modes), means, cis, 'k','linestyle','none','linewidth',1.2);
set(gca,'XTickLabel',modes,'FontName','Times New Roman','FontSize',12);
ylabel('Terminal true-vs-reference error (m)');
title('MPC endpoint accuracy by feedback source'); grid on
saveas(f1,'mc_terminal_error.png');

% (2) example trajectories (trial 1), true paths for each mode
f2 = figure('Color','w','Name','Example trajectories');
colors = {'k','r','b'}; hold on
plot(cfg.lm(1,:),cfg.lm(2,:),'k^','MarkerFaceColor','y','MarkerSize',9,'DisplayName','landmarks');
plot(cfg.xs(1),cfg.xs(2),'gp','MarkerFaceColor','g','MarkerSize',14,'DisplayName','goal');
for mi = 1:numel(modes)
    tr = results.(modes{mi}).sample.traj_true;
    plot(tr(1,:),tr(2,:),colors{mi},'linewidth',1.6,'DisplayName',modes{mi});
end
axis equal; grid on; legend('Location','best');
set(gca,'FontName','Times New Roman','FontSize',12);
xlabel('x (m)'); ylabel('y (m)'); title('True trajectories (trial 1)');
saveas(f2,'mc_trajectories.png');

fprintf('Saved mc_terminal_error.png and mc_trajectories.png\n');
