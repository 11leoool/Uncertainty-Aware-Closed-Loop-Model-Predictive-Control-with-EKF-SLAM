% RUN_MONTECARLO_UM  Stage 1 of the UNKNOWN-MAP refactor: free-space
% Monte-Carlo comparison (oracle / odom / slam) under the start-frame
% convention with fully unknown landmarks.
%
% Differences from the surveyed-landmark study (non-obstacle case\):
%   * No initial pose offset: the start pose defines the frame (zero initial
%     pose uncertainty). The old random offset would be unobservable here.
%   * Landmarks unknown: first-observation initialisation, uninformative prior
%     (std 10 m), known data association ("featured" landmarks).
%   * Actuation noise raised to std 0.1 (var 0.01): drift is what
%     distinguishes odometry from SLAM in this framing; at the old 1e-6 both
%     would be near-perfect and the comparison vacuous.
%   * New metric: final landmark-map error (slam only).

clear; clc; close all;
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

% ---------------- configuration ----------------
cfg.lm   = [-0.5 1 -1.5;
             0.5 1  0.0];
cfg.xs   = [1.5; 1.5; 0];
cfg.x0_nom = [0; 0; 0];

cfg.lm_uninf_std = 10;              % uninformative landmark prior (std, m)
cfg.var_v = 0.01;  cfg.var_w = 0.01;
cfg.var_d = 0.01;  cfg.var_a = 0.01;

cfg.sim_tim = 20;
cfg.tol     = 0.05;

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
    TR = zeros(M,1); RR = zeros(M,1); LE = zeros(M,1); ST = zeros(M,1); ME = nan(M,1);
    sample_traj = [];
    fprintf('Running mode: %-7s ', mode);
    for k = 1:M
        res = mc_run_trial_um(mode, mpc, cfg, noises{k});
        TR(k) = res.term_true_ref;
        RR(k) = res.rmse_true_ref;
        LE(k) = res.mean_true_est;
        ST(k) = res.steps;
        ME(k) = res.map_err;
        if k == 1, sample_traj = res; end
    end
    results.(mode).term_true_ref = TR;
    results.(mode).rmse_true_ref = RR;
    results.(mode).loc_err       = LE;
    results.(mode).steps         = ST;
    results.(mode).map_err       = ME;
    results.(mode).sample        = sample_traj;
    fprintf('done.\n');
end

% ---------------- report ----------------
ci = @(v) 1.96*std(v)/sqrt(numel(v));
lines = {};
lines{end+1} = sprintf('=== Unknown-map free-space Monte-Carlo (M = %d) ===', M);
lines{end+1} = sprintf('%-8s | %-22s | %-22s | %-8s | %-10s', 'mode', ...
        'term true-vs-ref [m]', 'mean loc err [m]', 'steps', 'map err [m]');
lines{end+1} = repmat('-',1,86);
for mi = 1:numel(modes)
    m  = modes{mi};
    TR = results.(m).term_true_ref;
    LE = results.(m).loc_err;
    lines{end+1} = sprintf('%-8s | %6.4f +/- %6.4f      | %6.4f +/- %6.4f      | %5.1f   | %8.4f', ...
            m, mean(TR), ci(TR), mean(LE), ci(LE), mean(results.(m).steps), ...
            mean(results.(m).map_err,'omitnan')); %#ok<SAGROW>
end
txt = strjoin(lines, newline);
disp(txt);
fid = fopen('mc_um_results.txt','w'); fprintf(fid,'%s\n',txt); fclose(fid);

save('mc_um_results.mat','results','cfg','M');
fprintf('Saved mc_um_results.mat / .txt\n');

% ---------------- figures ----------------
f1 = figure('Color','w','Name','Terminal true-vs-ref error');
means = cellfun(@(m) mean(results.(m).term_true_ref), modes);
cis   = cellfun(@(m) ci(results.(m).term_true_ref),   modes);
bar(means,0.6); hold on
errorbar(1:numel(modes), means, cis, 'k','linestyle','none','linewidth',1.2);
set(gca,'XTickLabel',modes,'FontName','Times New Roman','FontSize',12);
ylabel('Terminal true-vs-reference error (m)');
title('Unknown-map framing: endpoint accuracy by feedback source'); grid on
saveas(f1,'mc_um_terminal_error.png');

f2 = figure('Color','w','Name','Example trajectories');
colors = {'k','r','b'}; hold on
plot(cfg.lm(1,:),cfg.lm(2,:),'k^','MarkerFaceColor','y','MarkerSize',9,'DisplayName','landmarks (true)');
plot(cfg.xs(1),cfg.xs(2),'gp','MarkerFaceColor','g','MarkerSize',14,'DisplayName','goal');
for mi = 1:numel(modes)
    tr = results.(modes{mi}).sample.traj_true;
    plot(tr(1,:),tr(2,:),colors{mi},'linewidth',1.6,'DisplayName',modes{mi});
end
axis equal; grid on; legend('Location','best');
set(gca,'FontName','Times New Roman','FontSize',12);
xlabel('x (m)'); ylabel('y (m)'); title('Unknown map: true trajectories (trial 1)');
saveas(f2,'mc_um_trajectories.png');

% (3) pose-uncertainty evolution (slam sample): grows then contracts as
% landmarks are initialised and re-observed
f3 = figure('Color','w','Name','Uncertainty evolution');
rs = results.slam.sample;
plot((1:numel(rs.sig_t))*mpc.T, rs.sig_t, 'b-', 'LineWidth',1.5);
grid on; xlabel('time (s)'); ylabel('\sigma_{max} (m)');
title('SLAM pose uncertainty (trial 1): unknown-map framing');
set(gca,'FontName','Times New Roman','FontSize',12);
saveas(f3,'mc_um_sigma.png');

fprintf('Saved mc_um_terminal_error.png, mc_um_trajectories.png, mc_um_sigma.png\n');
