% RUN_MONTECARLO_OBS  Obstacle-avoidance Monte-Carlo comparison of MPC feedback
% sources (oracle / odom / slam) with a disk obstacle on the path. The MPC plans
% avoidance from its FEEDBACK pose, but collisions are judged on the TRUE pose,
% so localization error has a SAFETY consequence.
%
% Stage A: deterministic keep-out margin (cfg.cov_aware = false).
% Stage B: set cfg.cov_aware = true (and cfg.gamma) to inflate the slam margin
%          by the pose covariance (chance constraint).

clear; clc; close all;
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

% ---------------- scenario ----------------
cfg.lm   = [-0.5 1 -1.5;        % surveyed landmarks
             0.5 1  0.0];
cfg.xs   = [1.5; 1.5; 0];
cfg.x0_nom = [0; 0; 0];

cfg.sigma_init_pos = 0.10;  cfg.sigma_init_th = 0.05;
cfg.lm_prior_std   = 0.10;
cfg.var_v = 1e-6;  cfg.var_w = 1e-6;
cfg.var_d = 0.01;  cfg.var_a = 0.01;
cfg.sim_tim = 30;  cfg.tol = 0.05;

cfg.safe_buffer = 0.06;         % fixed planning buffer (all modes) -> oracle collision-free
cfg.cov_aware = true;           % STAGE B: covariance-aware (chance-constraint) slam margin
cfg.gamma     = 2.0;            % chance-constraint factor: margin += gamma*sqrt(lambda_max(Sigma_xy))

% obstacle on the start->goal diagonal: [ox oy r_obs]
obstacle = [0.5 0.5 0.15];

% ---------------- controller ----------------
opt.N = 14; opt.T = 0.2; opt.rob_diam = 0.3;
opt.v_max = 0.6; opt.omega_max = pi/4;
opt.obs = obstacle; opt.xy_min = -2; opt.xy_max = 2;
fprintf('Building obstacle NMPC ...\n');
mpc = mc_build_mpc_obs(opt);

L = size(cfg.lm,2);
maxiter = round(cfg.sim_tim/mpc.T);

% ---------------- pre-draw shared noise ----------------
M = 50; modes = {'oracle','odom','slam'};
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

% ---------------- run ----------------
results = struct();
for mi = 1:numel(modes)
    mode = modes{mi};
    TR = zeros(M,1); COL = false(M,1); MC = zeros(M,1); ST = zeros(M,1);
    sample = [];
    fprintf('Running mode: %-7s ', mode);
    for k = 1:M
        res = mc_run_trial_obs(mode, mpc, cfg, noises{k});
        TR(k)  = res.term_true_ref;
        COL(k) = res.collided;
        MC(k)  = res.min_clear;
        ST(k)  = res.steps;
        if k == 1, sample = res; end
    end
    results.(mode).term_true_ref = TR;
    results.(mode).collided      = COL;
    results.(mode).min_clear     = MC;
    results.(mode).steps         = ST;
    results.(mode).sample        = sample;
    fprintf('done.\n');
end

% ---------------- report ----------------
ci = @(v) 1.96*std(v)/sqrt(numel(v));
fprintf('\n=== Obstacle Monte-Carlo (M = %d, cov_aware = %d) ===\n', M, cfg.cov_aware);
fprintf('%-8s | %-20s | %-12s | %-18s\n','mode','term true-vs-ref [m]','collision %','min clearance [m]');
fprintf('%s\n', repmat('-',1,70));
for mi = 1:numel(modes)
    m  = modes{mi};
    TR = results.(m).term_true_ref;
    fprintf('%-8s | %6.4f +/- %6.4f    | %6.1f      | %+7.4f (mean %+.4f)\n', ...
        m, mean(TR), ci(TR), 100*mean(results.(m).collided), ...
        min(results.(m).min_clear), mean(results.(m).min_clear));
end
fprintf('%s\n', repmat('-',1,70));

save('mc_results_obs_stageB.mat','results','cfg','obstacle','M');
fprintf('Saved mc_results_obs_stageB.mat\n');

% ---------------- figures ----------------
ang = linspace(0,2*pi,100);

% (1) example avoidance trajectories (trial 1)
f1 = figure('Color','w','Name','Obstacle avoidance trajectories'); hold on
% obstacle disk + keep-out ring
ox=obstacle(1); oy=obstacle(2); orad=obstacle(3); safe=orad+mpc.rob_r;
fill(ox+orad*cos(ang), oy+orad*sin(ang), [0.85 0.5 0.5], 'EdgeColor','none');
plot(ox+safe*cos(ang), oy+safe*sin(ang), '--', 'Color',[0.6 0.2 0.2]);
plot(cfg.lm(1,:),cfg.lm(2,:),'k^','MarkerFaceColor','y','MarkerSize',9);
plot(cfg.xs(1),cfg.xs(2),'gp','MarkerFaceColor','g','MarkerSize',14);
colors = {'k','r','b'};
for mi = 1:numel(modes)
    tr = results.(modes{mi}).sample.traj_true;
    plot(tr(1,:),tr(2,:),colors{mi},'linewidth',1.6,'DisplayName',modes{mi});
end
legend({'obstacle','keep-out','landmarks','goal','oracle','odom','slam'},'Location','northwest');
axis equal; grid on; xlabel('x (m)'); ylabel('y (m)');
title('Stage B (covariance-aware): true trajectories (trial 1)');
set(gca,'FontName','Times New Roman','FontSize',12);
saveas(f1,'mc_obs_stageB_trajectories.png');

% (2) collision rate by mode
f2 = figure('Color','w','Name','Collision rate');
rates = cellfun(@(m) 100*mean(results.(m).collided), modes);
bar(rates,0.6);
set(gca,'XTickLabel',modes,'FontName','Times New Roman','FontSize',12);
ylabel('Collision rate (%)'); title('Stage B: collisions by feedback source'); grid on
saveas(f2,'mc_obs_stageB_collision.png');

fprintf('Saved mc_obs_stageB_trajectories.png and mc_obs_stageB_collision.png\n');
