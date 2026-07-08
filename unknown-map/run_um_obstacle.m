% RUN_UM_OBSTACLE  Stage 2 of the unknown-map refactor: waypoint tour with a
% STATIC obstacle, both map and obstacle UNKNOWN UNTIL DETECTED (finite
% sensing range r_sense). The obstacle sits on leg 2 of the tour, so it is
% encountered after drift has had a chance to accumulate.
%
% Modes (shared seeds rng(2024)):
%   oracle      true pose + clairvoyant true obstacle (upper bound)
%   odom_fixed  dead reckoning + sensed obstacle + fixed margin d0
%   slam_fixed  EKF-SLAM      + sensed obstacle + fixed margin d0
%   slam_cov    EKF-SLAM      + sensed obstacle + covariance-aware margin
%               delta = d0 + gamma*sqrt(lambda_max(P_xy + P_obs))

clear; clc; close all;
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

% ---------------- scenario ----------------
cfg.lm   = [-0.5 1 -1.5;
             0.5 1  0.0];
cfg.wps  = [ 1.5  -1.5   0.0;       % return-to-start patrol (matches Stage 1d)
             1.5   1.0   0.0];
cfg.x0_nom = [0; 0; 0];
cfg.obs  = [0 1.25 0.15];           % static obstacle on leg 2: [ox oy r]

cfg.lm_uninf_std = 10;
cfg.var_v = 0.01;  cfg.var_w = 0.01;
cfg.var_d = 0.01;  cfg.var_a = 0.01;

cfg.sim_tim = 45;
cfg.tol     = 0.05;
cfg.tol_wp  = 0.10;
cfg.wp_timeout = 15;                % s per intermediate waypoint (mission executive)
cfg.r_sense = 1.2;                  % landmarks AND obstacle unknown until detected

cfg.safe_buffer = 0.06;
cfg.gamma = 2.0;

% ---------------- controller (parameterised-obstacle NMPC) ----------------
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

% ---------------- run ----------------
modes = {'oracle','odom_fixed','slam_fixed','slam_cov'};
results = struct();
for mi = 1:numel(modes)
    mode = modes{mi};
    TR = zeros(M,1); COL = false(M,1); MC = zeros(M,1); PL = zeros(M,1);
    WP = zeros(M,1); DS = nan(M,1); OE = nan(M,1);
    sample = cell(M,1);
    fprintf('Running mode: %-11s ', mode);
    for k = 1:M
        res = mc_run_trial_um_obs(mode, mpc, cfg, noises{k});
        TR(k) = res.term_true_ref;  COL(k) = res.collided;
        MC(k) = res.min_clear;      PL(k)  = res.path_len;
        WP(k) = res.wp_reached;     DS(k)  = res.detect_step;
        OE(k) = res.obs_err_final;
        sample{k} = res;
    end
    results.(mode) = struct('term_true_ref',TR,'collided',COL,'min_clear',MC, ...
        'path_len',PL,'wp_reached',WP,'detect_step',DS,'obs_err_final',OE);
    results.(mode).sample = sample;
    fprintf('done (collisions: %d/%d).\n', sum(COL), M);
end

% ---------------- report ----------------
lines = {};
lines{end+1} = sprintf('=== Unknown-map obstacle study (M = %d, r_sense = %.1f, gamma = %.1f, d0 = %.2f) ===', ...
    M, cfg.r_sense, cfg.gamma, cfg.safe_buffer);
lines{end+1} = sprintf('%-11s | %-24s | %-10s | %-9s | %-9s | %-9s | %-8s', ...
    'mode','collision % (Wilson CI)','min clear','term err','path len','obs err','wps done');
lines{end+1} = repmat('-',1,104);
for mi = 1:numel(modes)
    m = modes{mi}; r = results.(m);
    nc_ = sum(r.collided);
    [lo,hi] = wilson(nc_, M);
    lines{end+1} = sprintf('%-11s | %5.1f%%  [%4.1f, %4.1f]%%  | %+8.4f  | %8.4f | %8.3f | %8.4f | %5.2f', ...
        m, 100*nc_/M, 100*lo, 100*hi, min(r.min_clear), mean(r.term_true_ref), ...
        mean(r.path_len), mean(r.obs_err_final,'omitnan'), mean(r.wp_reached)); %#ok<SAGROW>
end
txt = strjoin(lines, newline);
disp(txt);
fid = fopen('um_obstacle_results.txt','w'); fprintf(fid,'%s\n',txt); fclose(fid);

save('um_obstacle_results.mat','results','cfg','M');
fprintf('\nSaved um_obstacle_results.mat / .txt\n');

% ---------------- figures ----------------
ang = linspace(0,2*pi,100);
ox = cfg.obs(1); oy = cfg.obs(2); orad = cfg.obs(3);

% (1) collision rate by mode
f1 = figure('Color','w');
rates = cellfun(@(m) 100*mean(results.(m).collided), modes);
bar(rates, 0.6);
set(gca,'XTickLabel',{'oracle','odom fixed','slam fixed','slam cov'}, ...
    'FontName','Times New Roman','FontSize',12);
ylabel('Collision rate (%)'); grid on
title('Unknown-map obstacle study: collisions by strategy');
saveas(f1,'um_obs_collision.png');

% (2) sample trajectories (first slam_fixed collision seed, else trial 1)
si = find(results.slam_fixed.collided, 1); if isempty(si), si = 1; end
f2 = figure('Color','w'); hold on
fill(ox+orad*cos(ang), oy+orad*sin(ang), [0.85 0.5 0.5], 'EdgeColor','none');
plot(cfg.lm(1,:),cfg.lm(2,:),'k^','MarkerFaceColor','y','MarkerSize',9);
plot(cfg.wps(1,:),cfg.wps(2,:),'gp','MarkerFaceColor','g','MarkerSize',12);
cols = {'k','r','m','b'};
for mi = 1:numel(modes)
    tr = results.(modes{mi}).sample{si}.traj_true;
    plot(tr(1,:),tr(2,:),cols{mi},'LineWidth',1.4);
end
legend({'obstacle','landmarks','waypoints','oracle','odom fixed','slam fixed','slam cov'}, ...
       'Location','eastoutside');
axis equal; grid on; xlabel('x (m)'); ylabel('y (m)');
title(sprintf('Unknown-map obstacle: true trajectories (trial %d)', si));
set(gca,'FontName','Times New Roman','FontSize',12);
saveas(f2,'um_obs_trajectories.png');

% (3) margin & uncertainties over time (slam_cov, same trial)
rc = results.slam_cov.sample{si};
f3 = figure('Color','w'); hold on
t_c = (1:numel(rc.sig_t))*mpc.T;
plot(t_c, rc.sig_t, 'b-', 'LineWidth',1.4);
plot(t_c, rc.delta_t, 'b--', 'LineWidth',1.4);
plot(t_c, rc.obserr_t, 'k:', 'LineWidth',1.4);
grid on; xlabel('time (s)'); ylabel('(m)');
legend({'\sigma_{max} robot','applied inflation','obstacle est. error'}, ...
       'Location','northwest');
title(sprintf('slam\\_cov margin dynamics (trial %d)', si));
set(gca,'FontName','Times New Roman','FontSize',12);
saveas(f3,'um_obs_margin.png');

fprintf('Saved um_obs_collision.png, um_obs_trajectories.png, um_obs_margin.png\n');

% =========================== local functions ===============================
function [lo,hi] = wilson(k, n)
z = 1.96; p = k/n;
den = 1 + z^2/n;
ctr = (p + z^2/(2*n)) / den;
hw  = z*sqrt(p*(1-p)/n + z^2/(4*n^2)) / den;
lo = max(0, ctr-hw); hi = min(1, ctr+hw);
end
