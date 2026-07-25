% RUN_UM_RANDGEOM  Randomized-geometry Monte Carlo (static obstacle).
% Tests whether gamma = 2 -- tuned on the ORIGINAL fixed geometry -- transfers
% to layouts it never saw. Per geometry the obstacle position and the three
% landmark positions are randomized under feasibility constraints, then an
% ORACLE FEASIBILITY FILTER keeps only geometries the ideal-state controller
% solves collision-free (so any collision in the other modes is an
% estimation/margin failure, not an impossible scenario). gamma is FROZEN.
%
% Modes per accepted geometry: oracle / odom_fixed / slam_fixed / slam_cov,
% each on M shared noise seeds. Pooled + per-geometry collision rates.

clear; clc; close all;
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

% ---------------- base scenario (patrol; gamma FROZEN at the tuned value) ----------------
base.wps  = [ 1.5  -1.5   0.0;  1.5   1.0   0.0];
base.x0_nom = [0; 0; 0];
base.lm_uninf_std = 10;
base.var_v = 0.01;  base.var_w = 0.01;
base.var_d = 0.01;  base.var_a = 0.01;
base.sim_tim = 45;  base.tol = 0.05;  base.tol_wp = 0.10;  base.wp_timeout = 15;
base.r_sense = 1.2; base.safe_buffer = 0.06;  base.gamma = 2.0;
r_obs = 0.15; contact = 0.15 + r_obs;                 % robot + obstacle radii
L = 3;

% arena + path (for sampling)
arena = [-2 2; -1.5 2];                                 % x-range ; y-range
path = [base.x0_nom(1:2), base.wps, base.x0_nom(1:2)];  % 2 x 4 waypoint chain

opt.N=14; opt.T=0.2; opt.rob_diam=0.3; opt.v_max=0.6; opt.omega_max=pi/4;
opt.xy_min=-2; opt.xy_max=2;
fprintf('Building NMPC ...\n'); mpc = mc_build_mpc_dyn(opt);
maxit = round(base.sim_tim/mpc.T);

N_GEOM = 30; M = 20; MAXATT = 400;

% ---------------- pre-draw M shared noise realisations ----------------
rng(2024);
noises = cell(M,1);
for k=1:M
    nz.u = [sqrt(base.var_v); sqrt(base.var_w)].*randn(2,maxit);
    nz.z = zeros(2,L,maxit);
    nz.z(1,:,:) = sqrt(base.var_d)*randn(1,L,maxit);
    nz.z(2,:,:) = sqrt(base.var_a)*randn(1,L,maxit);
    nz.zo = [sqrt(base.var_d); sqrt(base.var_a)].*randn(2,maxit);
    noises{k}=nz;
end

% ---------------- sample + feasibility-filter geometries ----------------
fprintf('Sampling geometries (target %d) ...\n', N_GEOM);
geom = {}; att = 0;
while numel(geom) < N_GEOM && att < MAXATT
    att = att + 1;
    [lm, obs] = sample_geom(path, arena, r_obs, base.r_sense);
    if isempty(lm), continue; end
    c = base; c.lm = lm; c.obs = obs;
    feasible = true;                                    % oracle must solve it
    for s = 1:3
        r = mc_run_trial_um_obs('oracle', mpc, c, noises{s});
        if r.collided || r.term_true_ref > 0.15, feasible = false; break; end
    end
    if feasible, geom{end+1} = struct('lm',lm,'obs',obs); end %#ok<SAGROW>
end
NG = numel(geom);
fprintf('Accepted %d / %d geometries (%d attempts).\n', NG, N_GEOM, att);

% ---------------- run all modes on accepted geometries ----------------
modes = {'oracle','odom_fixed','slam_fixed','slam_cov'};
COLL = zeros(NG, numel(modes));       % collisions per geometry per mode (out of M)
MINCL= inf(NG, numel(modes));         % worst clearance per geometry per mode
for g = 1:NG
    c = base; c.lm = geom{g}.lm; c.obs = geom{g}.obs;
    fprintf('geom %2d/%d ', g, NG);
    for mi = 1:numel(modes)
        nc = 0; mcw = inf;
        for k = 1:M
            r = mc_run_trial_um_obs(modes{mi}, mpc, c, noises{k});
            nc = nc + r.collided; mcw = min(mcw, r.min_clear);
        end
        COLL(g,mi) = nc; MINCL(g,mi) = mcw;
    end
    fprintf('coll [orc odf slf slc] = %d %d %d %d\n', COLL(g,:));
end

% ---------------- aggregate ----------------
tot = NG*M;
lines = {};
lines{end+1} = sprintf('=== Randomized-geometry MC (static obstacle, gamma FROZEN=2) ===', 1);
lines{end+1} = sprintf('%d feasible geometries x %d seeds = %d trials/mode', NG, M, tot);
lines{end+1} = '';
lines{end+1} = sprintf('%-12s | %-22s | %-12s | %-14s', 'mode','pooled collision %% (CI)','worst clear','geoms w/ >=1 coll');
lines{end+1} = repmat('-',1,70);
for mi = 1:numel(modes)
    nc = sum(COLL(:,mi)); [lo,hi] = wil(nc, tot);
    ngc = sum(COLL(:,mi) > 0);
    lines{end+1} = sprintf('%-12s | %5.2f%% [%.2f,%.2f]      | %+8.4f    | %d / %d', ...
        modes{mi}, 100*nc/tot, 100*lo, 100*hi, min(MINCL(:,mi)), ngc, NG); %#ok<SAGROW>
end
txt = strjoin(lines, newline);
disp(txt);
fid=fopen('um_randgeom_results.txt','w'); fprintf(fid,'%s\n',txt); fclose(fid);
save('um_randgeom_results.mat','geom','COLL','MINCL','modes','NG','M','base','-v7');
fprintf('\nSaved um_randgeom_results.txt / .mat\n');

% =========================== local functions ===============================
function [lm, obs] = sample_geom(path, arena, r_obs, r_sense)
% obstacle: near a random path segment (offset < contact so avoidance needed)
nseg = size(path,2)-1;
seg = randi(nseg); t = 0.2 + 0.6*rand;
p = path(:,seg)*(1-t) + path(:,seg+1)*t;
ang = 2*pi*rand; off = 0.28*rand;               % within contact radius of path
oc = p + off*[cos(ang); sin(ang)];
oc(1) = min(max(oc(1),arena(1,1)+0.2), arena(1,2)-0.2);
oc(2) = min(max(oc(2),arena(2,1)+0.2), arena(2,2)-0.2);
obs = [oc(1) oc(2) r_obs];
% landmarks: 3, >=0.45 from obstacle, pairwise >=0.5, >=1 within r_sense of start
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
    if ok && ~any(vecnorm(cand - [0;0]) <= r_sense), ok=false; end   % coverage at start
    if ok, lm = cand; return; end
end
lm = [];   % failed -> caller skips
end

function [lo,hi]=wil(k,n)
z=1.96;p=k/n;den=1+z^2/n;ctr=(p+z^2/(2*n))/den;
hw=z*sqrt(p*(1-p)/n+z^2/(4*n^2))/den;lo=max(0,ctr-hw);hi=min(1,ctr+hw);
end
