% RUN_MONTECARLO_UM4  Stage 1d: RETURN-TO-START PATROL tour + sensing range 1.2 m. Stage 1c (r_sense 1.5, one-way tour) showed the final goal (1.2,-0.8) is landmark-denied (all lms >= 1.8 m) -> terminal metric degraded by construction. The patrol closes the loop back to (0,0), inside lm1 coverage: the canonical SLAM demonstration. Based on of the UNKNOWN-MAP refactor: multi-waypoint
% TOUR (~8.4 m) instead of the single 2.1 m dash of run_montecarlo_um.
%
% Rationale: in the start-frame convention there is no initial offset for
% SLAM to correct (Stage 1a showed slam 0.064 vs odom 0.071 m -- overlapping
% CIs); what SLAM buys is DRIFT SUPPRESSION, which is only material over
% distance. The tour makes dead-reckoning drift compound across three legs
% while SLAM keeps re-anchoring on the landmarks it mapped in leg 1.

clear; clc; close all;
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

% ---------------- configuration ----------------
cfg.lm   = [-0.5 1 -1.5;
             0.5 1  0.0];
cfg.wps  = [ 1.5  -1.5   0.0;          % waypoint x (visited in order)
             1.5   1.0   0.0];         % waypoint y ; last = return to start
cfg.x0_nom = [0; 0; 0];

cfg.lm_uninf_std = 10;
cfg.var_v = 0.01;  cfg.var_w = 0.01;
cfg.var_d = 0.01;  cfg.var_a = 0.01;

cfg.sim_tim = 45;
cfg.tol     = 0.05;             % final-goal arrival tolerance (terminal metric)
cfg.tol_wp  = 0.10;             % intermediate waypoints are pass-through
cfg.wp_timeout = 15;            % s per intermediate waypoint (mission executive)
cfg.r_sense = 1.2;              % landmark sensing range (m): unknown until detected

M     = 50;
modes = {'oracle','odom','slam'};

% ---------------- build controller once ----------------
fprintf('Building NMPC ...\n');
mpc = mc_build_mpc();
L = size(cfg.lm,2);
maxiter = round(cfg.sim_tim/mpc.T);

% ---------------- pre-draw noise (shared across modes) ----------------
rng(2024);
noises = cell(M,1);
for k = 1:M
    nz.u = [sqrt(cfg.var_v); sqrt(cfg.var_w)] .* randn(2, maxiter);
    nz.z = zeros(2, L, maxiter);
    nz.z(1,:,:) = sqrt(cfg.var_d) * randn(1, L, maxiter);
    nz.z(2,:,:) = sqrt(cfg.var_a) * randn(1, L, maxiter);
    noises{k} = nz;
end

% ---------------- run ----------------
results = struct();
for mi = 1:numel(modes)
    mode = modes{mi};
    TR = zeros(M,1); RR = zeros(M,1); LE = zeros(M,1); ST = zeros(M,1);
    ME = nan(M,1);  WP = zeros(M,1);
    sample_traj = [];
    fprintf('Running mode: %-7s ', mode);
    for k = 1:M
        res = mc_run_trial_um2(mode, mpc, cfg, noises{k});
        TR(k) = res.term_true_ref;
        RR(k) = res.rmse_true_ref;
        LE(k) = res.mean_true_est;
        ST(k) = res.steps;
        ME(k) = res.map_err;
        WP(k) = res.wp_reached;
        if k == 1, sample_traj = res; end
    end
    results.(mode).term_true_ref = TR;
    results.(mode).rmse_true_ref = RR;
    results.(mode).loc_err       = LE;
    results.(mode).steps         = ST;
    results.(mode).map_err       = ME;
    results.(mode).wp_reached    = WP;
    results.(mode).sample        = sample_traj;
    fprintf('done.\n');
end

% ---------------- report ----------------
ci = @(v) 1.96*std(v)/sqrt(numel(v));
lines = {};
lines{end+1} = sprintf('=== Unknown-map PATROL + sensing range (M = %d, r_sense = 1.2 m) ===', M);
lines{end+1} = sprintf('%-8s | %-22s | %-22s | %-8s | %-10s | %-8s', 'mode', ...
        'term true-vs-ref [m]', 'mean loc err [m]', 'steps', 'map err [m]', 'wps done');
lines{end+1} = repmat('-',1,98);
for mi = 1:numel(modes)
    m  = modes{mi};
    TR = results.(m).term_true_ref;
    LE = results.(m).loc_err;
    lines{end+1} = sprintf('%-8s | %6.4f +/- %6.4f      | %6.4f +/- %6.4f      | %5.1f   | %8.4f   | %4.2f', ...
            m, mean(TR), ci(TR), mean(LE), ci(LE), mean(results.(m).steps), ...
            mean(results.(m).map_err,'omitnan'), mean(results.(m).wp_reached)); %#ok<SAGROW>
end
txt = strjoin(lines, newline);
disp(txt);
fid = fopen('mc_um4_results.txt','w'); fprintf(fid,'%s\n',txt); fclose(fid);

save('mc_um4_results.mat','results','cfg','M');
fprintf('Saved mc_um4_results.mat / .txt\n');

% ---------------- figures ----------------
f1 = figure('Color','w');
means = cellfun(@(m) mean(results.(m).term_true_ref), modes);
cis   = cellfun(@(m) ci(results.(m).term_true_ref),   modes);
bar(means,0.6); hold on
errorbar(1:numel(modes), means, cis, 'k','linestyle','none','linewidth',1.2);
set(gca,'XTickLabel',modes,'FontName','Times New Roman','FontSize',12);
ylabel('Terminal true-vs-reference error (m)');
title('Unknown-map tour: endpoint accuracy by feedback source'); grid on
saveas(f1,'mc_um4_terminal_error.png');

f2 = figure('Color','w');
colors = {'k','r','b'}; hold on
plot(cfg.lm(1,:),cfg.lm(2,:),'k^','MarkerFaceColor','y','MarkerSize',9,'DisplayName','landmarks (true)');
plot(cfg.wps(1,:),cfg.wps(2,:),'gp','MarkerFaceColor','g','MarkerSize',12,'DisplayName','waypoints');
for mi = 1:numel(modes)
    tr = results.(modes{mi}).sample.traj_true;
    plot(tr(1,:),tr(2,:),colors{mi},'linewidth',1.4,'DisplayName',modes{mi});
end
axis equal; grid on; legend('Location','best');
set(gca,'FontName','Times New Roman','FontSize',12);
xlabel('x (m)'); ylabel('y (m)'); title('Unknown-map tour: true trajectories (trial 1)');
saveas(f2,'mc_um4_trajectories.png');

f3 = figure('Color','w'); hold on
rs = results.slam.sample;
ro = results.odom.sample;
plot((1:numel(ro.sig_t))*mpc.T, ro.sig_t, 'r-.', 'LineWidth',1.4);
plot((1:numel(rs.sig_t))*mpc.T, rs.sig_t, 'b-',  'LineWidth',1.4);
grid on; xlabel('time (s)'); ylabel('\sigma_{max} (m)');
legend({'odom (dead reckoning)','slam'},'Location','northwest');
title('Pose uncertainty over the tour (trial 1)');
set(gca,'FontName','Times New Roman','FontSize',12);
saveas(f3,'mc_um4_sigma.png');

fprintf('Saved mc_um4_terminal_error.png, mc_um4_trajectories.png, mc_um4_sigma.png\n');



