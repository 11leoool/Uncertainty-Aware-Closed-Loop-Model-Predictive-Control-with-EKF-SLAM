% RUN_UM_CONSISTENCY  Filter-consistency stress test (robustness study).
%
% The covariance-aware margin trusts the filter's reported covariance. This
% sweep asks how much the filter may be MIS-CALIBRATED before safety breaks.
% Knob: filt_scale = (filter's assumed noise variance) / (true noise variance).
% The TRUE noise draws are held fixed (shared seeds); only the filter's belief
% changes, so the point estimate and reported covariance shift together as a
% genuinely (in)consistent filter would.
%   filt_scale < 1 : OVERCONFIDENT (reports P smaller than the true error)
%   filt_scale = 1 : consistent (the paper's primary setting)
%   filt_scale > 1 : conservative (over-reports uncertainty)
%
% Reported per scale: cov-aware collision rate (Wilson CI), min clearance,
% mean reported sigma vs mean TRUE localization error (the consistency gap),
% path. slam_fixed (delta0, covariance-blind) shown once as a reference floor.
% Scenario = sensed static obstacle (run_um_obstacle.m configuration).

clear; clc; close all;
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

cfg.lm   = [-0.5 1 -1.5;
             0.5 1  0.0];
cfg.wps  = [ 1.5  -1.5   0.0;
             1.5   1.0   0.0];
cfg.x0_nom = [0; 0; 0];
cfg.obs  = [0 1.25 0.15];
cfg.lm_uninf_std = 10;
cfg.var_v = 0.01;  cfg.var_w = 0.01;
cfg.var_d = 0.01;  cfg.var_a = 0.01;
cfg.sim_tim = 45;  cfg.tol = 0.05;  cfg.tol_wp = 0.10;  cfg.wp_timeout = 15;
cfg.r_sense = 1.2; cfg.safe_buffer = 0.06;  cfg.gamma = 2.0;

SCALES = [0.1 0.25 0.5 1.0 2.0 4.0];   % filter assumed / true noise variance

opt.N = 14; opt.T = 0.2; opt.rob_diam = 0.3;
opt.v_max = 0.6; opt.omega_max = pi/4;
opt.xy_min = -2; opt.xy_max = 2;
fprintf('Building parameterised-obstacle NMPC ...\n');
mpc = mc_build_mpc_dyn(opt);

L = size(cfg.lm,2);
maxiter = round(cfg.sim_tim/mpc.T);

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

% reference: covariance-blind fixed margin (filt_scale irrelevant, consistent)
ref = run_cell('slam_fixed', 1.0, mpc, cfg, noises, M);

nS = numel(SCALES);
sweep = struct('scale',SCALES,'coll',zeros(1,nS),'lo',zeros(1,nS),'hi',zeros(1,nS), ...
    'minclear',zeros(1,nS),'sigma',zeros(1,nS),'locerr',zeros(1,nS),'path',zeros(1,nS));
lines = {};
lines{end+1} = sprintf('=== Filter-consistency sweep (M = %d, static sensed obstacle, cov-aware gamma=2) ===', M);
lines{end+1} = sprintf('reference slam_fixed (delta0, cov-blind): %s', ref.str);
lines{end+1} = '';
lines{end+1} = sprintf('%-8s | %-22s | %-10s | %-10s | %-10s | %-8s', ...
    'filt/true','cov coll %% (Wilson)','min clear','mean sigma','TRUE locerr','path');
lines{end+1} = repmat('-',1,86);
for i = 1:nS
    c = run_cell('slam_cov', SCALES(i), mpc, cfg, noises, M);
    sweep.coll(i)=c.coll; sweep.lo(i)=c.lo; sweep.hi(i)=c.hi;
    sweep.minclear(i)=c.minclear; sweep.sigma(i)=c.sigma;
    sweep.locerr(i)=c.locerr; sweep.path(i)=c.path;
    tag = ''; if SCALES(i)<1, tag=' (overconf)'; elseif SCALES(i)>1, tag=' (conserv)'; else, tag=' (consistent)'; end
    lines{end+1} = sprintf('%6.2f%-9s| %5.1f%% [%4.1f,%4.1f]      | %+8.4f  | %8.4f   | %8.4f    | %6.3f', ...
        SCALES(i), tag, c.coll, c.lo, c.hi, c.minclear, c.sigma, c.locerr, c.path); %#ok<SAGROW>
end
txt = strjoin(lines, newline);
disp(txt);
fid = fopen('um_consistency_results.txt','w'); fprintf(fid,'%s\n',txt); fclose(fid);
save('um_consistency_results.mat','sweep','ref','SCALES','cfg','M');

% ---- figure: collision rate + consistency gap vs filt_scale (log x) ----
f = figure('Color','w');
yyaxis left
errorbar(SCALES, sweep.coll, sweep.coll-sweep.lo, sweep.hi-sweep.coll, ...
    '-o','LineWidth',1.5,'MarkerFaceColor','auto'); hold on
plot(SCALES, ref.coll*ones(size(SCALES)), '--', 'LineWidth',1.0);
set(gca,'XScale','log'); ylabel('Collision rate (%)'); ylim([-2 max(sweep.coll)*1.15+2]);
yyaxis right
plot(SCALES, sweep.sigma, '-s', 'LineWidth',1.3);
plot(SCALES, sweep.locerr, ':^', 'LineWidth',1.3);
ylabel('reported \sigma / true loc. error (m)');
xlabel('filter-assumed / true noise variance');
legend({'cov-aware collisions','fixed-margin ref','reported \sigma','true loc. error'}, ...
    'Location','north','FontSize',7);
grid on; title('Robustness to filter (in)consistency');
set(gca,'FontName','Times New Roman','FontSize',11);
saveas(f,'um_consistency.png');
fprintf('\nSaved um_consistency_results.txt/.mat and um_consistency.png\n');

% =========================== local functions ===============================
function r = run_cell(md, fs, mpc, cfg, noises, M)
c = cfg; c.filt_scale = fs;
COL=false(M,1); MC=zeros(M,1); SG=zeros(M,1); LE=zeros(M,1); PL=zeros(M,1);
fprintf('Running %-10s filt_scale=%.2f ', md, fs);
for k = 1:M
    res = mc_run_trial_um_obs(md, mpc, c, noises{k});
    COL(k)=res.collided; MC(k)=res.min_clear; SG(k)=res.mean_sigma;
    LE(k)=res.mean_loc_err; PL(k)=res.path_len;
end
[lo,hi] = wilson(sum(COL), M);
r.coll=100*mean(COL); r.lo=100*lo; r.hi=100*hi; r.minclear=min(MC);
r.sigma=mean(SG); r.locerr=mean(LE); r.path=mean(PL);
r.str = sprintf('coll %.1f%% [%.1f,%.1f], min clear %+.4f, path %.3f', ...
    r.coll, r.lo, r.hi, r.minclear, r.path);
fprintf('-> coll %.1f%%\n', r.coll);
end

function [lo,hi] = wilson(k, n)
z = 1.96; p = k/n;
den = 1 + z^2/n;
ctr = (p + z^2/(2*n)) / den;
hw  = z*sqrt(p*(1-p)/n + z^2/(4*n^2)) / den;
lo = max(0, ctr-hw); hi = min(1, ctr+hw);
end
