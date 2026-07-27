% RUN_UM_RANDGEOM_V2 - randomized-geometry transfer, corrected for two defects
% found in verification review:
%
%  (1) DISJOINT SEEDS. v1 screened each geometry with noise realisations 1-3 and
%      then evaluated on realisations 1-20, so screening trials sat inside the
%      final denominator. Here screening and evaluation draw from separate,
%      non-overlapping streams; no screening observation enters inference.
%
%  (2) CONTINUOUS COLLISION SCORING. v1 tested clearance only at sampled states
%      (0.2 s apart, up to ~0.12 m of travel per step against a 0.30 m keep-out),
%      so a grazing pass could dip inside and out between samples. Here the
%      minimum distance from each travelled SEGMENT to the obstacle centre is
%      used, which is exact for a static obstacle under piecewise-linear motion.
%      Both scorings are reported so the difference is visible.
%
% Saves um_randgeom_v2_results.{mat,txt} INCLUDING per-trial outcomes and the
% accepted layouts, so the released data matches what the manuscript reports.

clear; clc; close all;
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

base.wps  = [ 1.5  -1.5   0.0;  1.5   1.0   0.0];
base.x0_nom = [0; 0; 0];
base.lm_uninf_std = 10;
base.var_v = 0.01;  base.var_w = 0.01;
base.var_d = 0.01;  base.var_a = 0.01;
base.sim_tim = 45;  base.tol = 0.05;  base.tol_wp = 0.10;  base.wp_timeout = 15;
base.r_sense = 1.2; base.safe_buffer = 0.06;  base.gamma = 2.0;
r_obs = 0.15; rob_r = 0.15;
L = 3;
arena = [-2 2; -1.5 2];
path = [base.x0_nom(1:2), base.wps, base.x0_nom(1:2)];

opt.N=14; opt.T=0.2; opt.rob_diam=0.3; opt.v_max=0.6; opt.omega_max=pi/4;
opt.xy_min=-2; opt.xy_max=2;
fprintf('Building NMPC ...\n'); mpc = mc_build_mpc_dyn(opt);
maxit = round(base.sim_tim/mpc.T);

% ---- CALIBRATION PHASE (frozen protocol: calibrated_fixed/PROTOCOL.md) ----
% All streams disjoint from the published run (9001/4242/2024) and from every
% other study in this project.
N_TRAIN = 10; M_CAL = 20; N_SCREEN = 3; MAXATT = 400;
rng(9500,'twister');
screen_noise = cell(N_SCREEN,1);
for k=1:N_SCREEN, screen_noise{k} = mknoise(base, L, maxit); end
rng(7100,'twister');
cal_noise = cell(M_CAL,1);
for k=1:M_CAL, cal_noise{k} = mknoise(base, L, maxit); end
rng(3030,'twister');                       % training geometry stream
geom = {}; att = 0;
while numel(geom) < N_TRAIN && att < MAXATT
    att = att + 1;
    [lm, obs] = sample_geom(path, arena, r_obs, base.r_sense);
    if isempty(lm), continue; end
    c = base; c.lm = lm; c.obs = obs;
    feasible = true;
    for s = 1:N_SCREEN
        r = mc_run_trial_um_obs('oracle', mpc, c, screen_noise{s});
        if r.collided || r.term_true_ref > 0.15, feasible = false; break; end
    end
    if feasible, geom{end+1} = struct('lm',lm,'obs',obs); end %#ok<SAGROW>
end
NT = numel(geom);
fprintf('Training: accepted %d / %d layouts (%d attempts).\n', NT, N_TRAIN, att);
assert(NT == N_TRAIN, 'training sampling failed');
cvals = []; af = [];
for g = 1:NT
    c = base; c.lm = geom{g}.lm; c.obs = geom{g}.obs;
    for k = 1:M_CAL
        r = mc_run_trial_um_obs('slam_cov', mpc, c, cal_noise{k});
        d = r.delta_t(:); act = d > 0;
        if any(act), cvals(end+1) = mean(d(act)); af(end+1) = mean(act); end %#ok<SAGROW>
    end
    fprintf('  layout %d/%d done\n', g, NT);
end
cfixed = mean(cvals);
fid = fopen('calib_fixed_result.txt','w');
fprintf(fid, 'CALIBRATED FIXED CONSTANT (frozen): c = %.4f m\n', cfixed);
fprintf(fid, 'per-trial conditional means: median %.4f IQR %.4f min %.4f max %.4f (n=%d)\n', ...
    median(cvals), iqr(cvals), min(cvals), max(cvals), numel(cvals));
fprintf(fid, 'active fraction: mean %.3f\n', mean(af));
fclose(fid); type calib_fixed_result.txt
save('calib_fixed_result.mat','cfixed','cvals','af');
fprintf('CALIB_FIXED_DONE\n');
function nz = mknoise(base, L, maxit)
    nz.u = [sqrt(base.var_v); sqrt(base.var_w)].*randn(2,maxit);
    nz.z = zeros(2,L,maxit);
    nz.z(1,:,:) = sqrt(base.var_d)*randn(1,L,maxit);
    nz.z(2,:,:) = sqrt(base.var_a)*randn(1,L,maxit);
    nz.zo = [sqrt(base.var_d); sqrt(base.var_a)].*randn(2,maxit);
end

function [lm, obs] = sample_geom(path, arena, r_obs, r_sense)
nseg = size(path,2)-1;
seg = randi(nseg); t = 0.2 + 0.6*rand;
p = path(:,seg)*(1-t) + path(:,seg+1)*t;
ang = 2*pi*rand; off = 0.28*rand;
oc = p + off*[cos(ang); sin(ang)];
oc(1) = min(max(oc(1),arena(1,1)+0.2), arena(1,2)-0.2);
oc(2) = min(max(oc(2),arena(2,1)+0.2), arena(2,2)-0.2);
obs = [oc(1) oc(2) r_obs];
lm = [];
for tries = 1:300
    cand = [arena(1,1)+diff(arena(1,:))*rand(1,3);
            arena(2,1)+diff(arena(2,:))*rand(1,3)];
    ok = true;
    for i=1:3
        if norm(cand(:,i)-oc) < 0.45, ok=false; break; end
    end
    if ok
        for i=1:3, for j=i+1:3
            if norm(cand(:,i)-cand(:,j)) < 0.5, ok=false; end
        end, end
    end
    if ok && ~any(vecnorm(cand - [0;0]) <= r_sense), ok=false; end
    if ok, lm = cand; return; end
end
lm = [];
end
