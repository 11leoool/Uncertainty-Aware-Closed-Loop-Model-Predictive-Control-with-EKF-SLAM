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

N_GEOM = 30; M_EVAL = 20; N_SCREEN = 3; MAXATT = 400;

% ---------- SEPARATE seed streams ----------
rng(9001,'twister');                       % screening stream (never evaluated)
screen_noise = cell(N_SCREEN,1);
for k=1:N_SCREEN, screen_noise{k} = mknoise(base, L, maxit); end
rng(4242,'twister');                       % evaluation stream (disjoint)
eval_noise = cell(M_EVAL,1);
for k=1:M_EVAL, eval_noise{k} = mknoise(base, L, maxit); end
rng(2024,'twister');                       % geometry sampling stream

fprintf('Sampling + screening geometries (screen seeds are NOT evaluated)...\n');
geom = {}; att = 0;
while numel(geom) < N_GEOM && att < MAXATT
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
NG = numel(geom);
fprintf('Accepted %d / %d geometries (%d attempts).\n', NG, N_GEOM, att);


% ---- EVALUATION of the frozen calibrated constant (PROTOCOL.md) ----
L0 = load('calibrated_fixed/calib_fixed_result.mat'); cfixed = L0.cfixed;
fprintf('frozen calibrated constant c = %.4f m\n', cfixed);
S = load('um_randgeom_v2_results.mat');

% replay validation: 3x3 cells of plain slam_fixed must match stored bit-exactly
fprintf('Replay validation...\n');
for g = 1:3
    c = base; c.lm = geom{g}.lm; c.obs = geom{g}.obs;
    for k = 1:3
        r = mc_run_trial_um_obs('slam_fixed', mpc, c, eval_noise{k});
        clc_ = cont_score(r, geom{g}.obs(1:2)', rob_r, r_obs);
        idx = find([S.TRIAL.geom]==g & [S.TRIAL.seed]==k & strcmp({S.TRIAL.mode},'slam_fixed'), 1);
        assert(~isempty(idx), 'stored record missing');
        assert(abs(clc_ - S.TRIAL(idx).min_clear_cont) < 1e-12, ...
            'REPLAY MISMATCH at g=%d k=%d: regenerated %g vs stored %g', ...
            g, k, clc_, S.TRIAL(idx).min_clear_cont);
    end
end
fprintf('Replay validation PASSED (9/9 cells bit-exact).\n');

R = struct('geom',{},'seed',{},'collc',{},'clr',{},'path',{});
tic;
for g = 1:NG
    c = base; c.lm = geom{g}.lm; c.obs = geom{g}.obs; c.fixed_extra = cfixed;
    for k = 1:M_EVAL
        r = mc_run_trial_um_obs('slam_fixed', mpc, c, eval_noise{k});
        clc_ = cont_score(r, geom{g}.obs(1:2)', rob_r, r_obs);
        R(end+1) = struct('geom',g,'seed',k,'collc',clc_<0, ...
                          'clr',clc_,'path',r.path_len); %#ok<SAGROW>
    end
    fprintf('  geom %d/%d (%.0fs)\n', g, NG, toc);
end
save('calibrated_fixed/eval_fixed_results.mat','R','cfixed');

% ---- paired two-way analysis vs stored slam_cov / slam_fixed ----
NGe = NG; MS = M_EVAL;
Ccal = zeros(NGe,MS); Ccov = zeros(NGe,MS); Cfix = zeros(NGe,MS);
Pcal = zeros(NGe,MS);
for i = 1:numel(R), Ccal(R(i).geom,R(i).seed) = R(i).collc; Pcal(R(i).geom,R(i).seed) = R(i).path; end  % stored TRIAL carries no path
for i = 1:numel(S.TRIAL)
    t = S.TRIAL(i);
    if strcmp(t.mode,'slam_cov'),   Ccov(t.geom,t.seed) = t.min_clear_cont<0; end
    if strcmp(t.mode,'slam_fixed'), Cfix(t.geom,t.seed) = t.min_clear_cont<0; end
end
fid = fopen('calibrated_fixed/eval_fixed_analysis.txt','w');
fprintf(fid, '=== Calibrated fixed constant (c = %.4f m, frozen on disjoint training) ===\n', cfixed);
fprintf(fid, 'collisions: cal-fixed %d/600 | slam_cov %d/600 | delta0-fixed %d/600\n', ...
    sum(Ccal(:)), sum(Ccov(:)), sum(Cfix(:)));
fprintf(fid, 'mean path:  cal-fixed %.3f (slam_cov path not stored per-trial in v2)\n', mean(Pcal(:)));
rng(11,'twister'); B = 20000;
D = Ccov - Ccal;                                   % negative = cov-aware better
bs = zeros(B,1);
for b = 1:B
    gi = randi(NGe,NGe,1); si = randi(MS,MS,1);
    sub = D(gi,si); bs(b) = mean(sub(:));
end
fprintf(fid, 'slam_cov - cal-fixed rate diff %+.4f  two-way 95%% CI [%+.4f,%+.4f]\n', ...
    mean(D(:)), quantile(bs,0.025), quantile(bs,0.975));
fprintf(fid, 'Pre-declared: comparable performance -> frame as automatic margin\n');
fprintf(fid, 'generation + transfer; small counts may be inconclusive and are\n');
fprintf(fid, 'reported as such (PROTOCOL.md).\n');
fclose(fid); type calibrated_fixed/eval_fixed_analysis.txt
fprintf('EVAL_FIXED_DONE\n');

function clc_ = cont_score(r, oc, rob_r, r_obs)
    P2 = r.traj_true(1:2,:); mind = inf;
    for j = 1:size(P2,2)-1
        pa = P2(:,j); pb = P2(:,j+1); dseg = pb - pa; L2 = dseg'*dseg;
        if L2 > 0, tt = max(0,min(1,((oc-pa)'*dseg)/L2)); else, tt = 0; end
        mind = min(mind, norm(pa + tt*dseg - oc));
    end
    clc_ = mind - (rob_r + r_obs);
end

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
