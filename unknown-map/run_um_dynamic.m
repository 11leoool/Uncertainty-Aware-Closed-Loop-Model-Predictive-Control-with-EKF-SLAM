% RUN_UM_DYNAMIC  Stage 4 of the unknown-map refactor: MOVING obstacle
% crossing patrol leg 1, unknown until detected (finite sensing range),
% tracked by a CV-EKF that coasts (prediction-only) when out of range.
% Four strategies on shared seeds; collisions vs the TRUE moving obstacle.

clear; clc; close all;
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

% ---------------- scenario (patrol as Stages 1d-3) ----------------
cfg.lm   = [-0.5 1 -1.5;
             0.5 1  0.0];
cfg.wps  = [ 1.5  -1.5   0.0;
             1.5   1.0   0.0];
cfg.x0_nom = [0; 0; 0];

cfg.lm_uninf_std = 10;
cfg.var_v = 0.01;  cfg.var_w = 0.01;
cfg.var_d = 0.01;  cfg.var_a = 0.01;

% Dynamic scenario uses a LONGER waypoint budget and mission time than the
% static studies: a safe controller must be allowed to wait out a crossing
% obstacle rather than being forced to abandon the waypoint (with 15 s the
% covariance-aware controller skipped wp1 in 45/50 trials). The timeout still
% bounds the hover-wedge pathology it was introduced for.
cfg.sim_tim = 75;
cfg.tol     = 0.05;
cfg.tol_wp  = 0.10;
cfg.wp_timeout = 30;
cfg.r_sense = 1.2;
cfg.safe_buffer = 0.06;
cfg.gamma = 2.0;

% moving obstacle: crosses leg 1 (the start->NE diagonal)
cfg.o0  = [1.3; 0.15];          % initial true position
cfg.vo0 = [-0.16; 0.16];        % initial true velocity (crossing)
cfg.obs_r = 0.15;
% True accel std 0.05: a moderately manoeuvring obstacle whose speed stays
% below the robot's over the 75 s mission. (The former 0.15, calibrated for
% 20-30 s missions, random-walked the obstacle to 0.46-0.85 m/s by mission
% end -- faster than the robot -- making safe waypoint visits impossible.)
% The filter's 0.02 over-estimates the true 0.0025: conservative coasting.
cfg.obs_sa2_true   = 0.0025;
cfg.obs_sa2_filter = 0.02;
cfg.track_drop = 3.0;           % s without detection -> drop stale track

% ---------------- controller ----------------
opt.N = 14; opt.T = 0.2; opt.rob_diam = 0.3;
opt.v_max = 0.6; opt.omega_max = pi/4;
opt.xy_min = -2; opt.xy_max = 2;
fprintf('Building parameterised-obstacle NMPC ...\n');
mpc = mc_build_mpc_dyn(opt);

L = size(cfg.lm,2);
maxiter = round(cfg.sim_tim/mpc.T);
ntot = maxiter + mpc.N + 2;

% ---------------- shared noise ----------------
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

% ---------------- run ----------------
modes = {'oracle','static','cv_fixed','cv_cov'};
results = struct();
for mi = 1:numel(modes)
    mode = modes{mi};
    TR = zeros(M,1); COL = false(M,1); MC = zeros(M,1); PL = zeros(M,1);
    WP = zeros(M,1); DS = nan(M,1); WV = zeros(M,1); LE = zeros(M,1);
    MD = zeros(M,1);
    sample = cell(M,1);
    fprintf('Running mode: %-9s ', mode);
    for k = 1:M
        res = mc_run_trial_um_dyn(mode, mpc, cfg, noises{k});
        TR(k) = res.term_true_ref;  COL(k) = res.collided;
        MC(k) = res.min_clear;      PL(k)  = res.path_len;
        WP(k) = res.wp_reached;     DS(k)  = res.detect_step;
        LE(k) = res.mean_loc_err;
        % waypoints TRULY visited (< 0.15 m) + mean true miss distance per wp
        wv = 0; md = zeros(1,size(cfg.wps,2));
        for w = 1:size(cfg.wps,2)
            md(w) = min(vecnorm(res.traj(1:2,:) - cfg.wps(:,w)));
            if md(w) < 0.15, wv = wv + 1; end
        end
        WV(k) = wv;  MD(k) = mean(md);
        sample{k} = res;
    end
    results.(mode) = struct('term_true_ref',TR,'collided',COL,'min_clear',MC, ...
        'path_len',PL,'wp_reached',WP,'detect_step',DS,'wp_visited',WV, ...
        'loc_err',LE,'wp_missdist',MD);
    results.(mode).sample = sample;
    fprintf('done (collisions: %d/%d, wps %.2f/3, locerr %.3f).\n', ...
        sum(COL), M, mean(WV), mean(LE));
end

% ---------------- report ----------------
lines = {};
lines{end+1} = sprintf('=== Unknown-map DYNAMIC obstacle (M = %d, r_sense = %.1f, gamma = %.1f) ===', ...
    M, cfg.r_sense, cfg.gamma);
lines{end+1} = sprintf('%-9s | %-24s | %-10s | %-9s | %-9s | %-11s | %-9s | %-9s', ...
    'mode','collision % (Wilson CI)','min clear','term err','path len','wps visited','loc err','wp missd');
lines{end+1} = repmat('-',1,112);
for mi = 1:numel(modes)
    m = modes{mi}; r = results.(m);
    nc_ = sum(r.collided);
    [lo,hi] = wilson(nc_, M);
    lines{end+1} = sprintf('%-9s | %5.1f%%  [%4.1f, %4.1f]%%  | %+8.4f  | %8.4f | %8.3f | %5.2f/3    | %8.4f | %8.4f', ...
        m, 100*nc_/M, 100*lo, 100*hi, min(r.min_clear), mean(r.term_true_ref), ...
        mean(r.path_len), mean(r.wp_visited), mean(r.loc_err), mean(r.wp_missdist)); %#ok<SAGROW>
end
txt = strjoin(lines, newline);
disp(txt);
fid = fopen('um_dynamic_results.txt','w'); fprintf(fid,'%s\n',txt); fclose(fid);
save('um_dynamic_results.mat','results','cfg','M');
fprintf('\nSaved um_dynamic_results.mat / .txt\n');

% ---------------- figures ----------------
f1 = figure('Color','w');
rates = cellfun(@(m) 100*mean(results.(m).collided), modes);
bar(rates, 0.6);
set(gca,'XTickLabel',{'oracle','static','cv fixed','cv cov'}, ...
    'FontName','Times New Roman','FontSize',12);
ylabel('Collision rate (%)'); grid on
title('Unknown-map dynamic obstacle: collisions by strategy');
saveas(f1,'um_dyn_collision.png');

si = find(results.static.collided, 1); if isempty(si), si = 1; end
f2 = figure('Color','w'); hold on
plot(cfg.lm(1,:),cfg.lm(2,:),'k^','MarkerFaceColor','y','MarkerSize',9);
plot(cfg.wps(1,:),cfg.wps(2,:),'gp','MarkerFaceColor','g','MarkerSize',12);
rs = results.static.sample{si};  rc = results.cv_cov.sample{si};
plot(rs.otraj(1,:), rs.otraj(2,:), 'Color',[0.85 0.5 0.5], 'LineWidth',2.5);
plot(rs.traj(1,:), rs.traj(2,:), 'r-.', 'LineWidth',1.5);
plot(rc.traj(1,:), rc.traj(2,:), 'b-',  'LineWidth',1.5);
legend({'landmarks','waypoints','obstacle (true path)','static (freeze)','cv\_cov'}, ...
       'Location','eastoutside');
axis equal; grid on; xlabel('x (m)'); ylabel('y (m)');
title(sprintf('Dynamic obstacle, unknown map (trial %d)', si));
set(gca,'FontName','Times New Roman','FontSize',12);
saveas(f2,'um_dyn_trajectories.png');

fprintf('Saved um_dyn_collision.png, um_dyn_trajectories.png\n');

% =========================== local functions ===============================
function [lo,hi] = wilson(k, n)
z = 1.96; p = k/n;
den = 1 + z^2/n;
ctr = (p + z^2/(2*n)) / den;
hw  = z*sqrt(p*(1-p)/n + z^2/(4*n^2)) / den;
lo = max(0, ctr-hw); hi = min(1, ctr+hw);
end
