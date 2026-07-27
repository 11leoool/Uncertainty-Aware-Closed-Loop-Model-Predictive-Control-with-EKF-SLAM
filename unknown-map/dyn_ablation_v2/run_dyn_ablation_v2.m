% RUN_UM_DYN_ABLATION  Dynamic matched-mean ablation (review roadmap item 1).
% Phase 1: run cv_cov (gamma = 2) and measure the mean node-0 inflation it
%          actually applies over all trials.
% Phase 2: run cv_fixed with a CONSTANT inflation equal to that mean, on
%          identical seeds. If the size-matched constant is equally safe, the
%          static ablation's conclusion (safety ~ margin size; adaptivity's
%          value = self-tuning) extends to the dynamic case; if not, the
%          adaptivity claim gains direct evidence. Either outcome is reported.
% Scenario identical to run_um_dynamic.m (final configuration).

clear; clc; close all;
addpath('..'); addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

cfg.lm   = [-0.5 1 -1.5;
             0.5 1  0.0];
cfg.wps  = [ 1.5  -1.5   0.0;
             1.5   1.0   0.0];
cfg.x0_nom = [0; 0; 0];
cfg.lm_uninf_std = 10;
cfg.var_v = 0.01;  cfg.var_w = 0.01;
cfg.var_d = 0.01;  cfg.var_a = 0.01;
cfg.sim_tim = 75;  cfg.tol = 0.05;  cfg.tol_wp = 0.10;  cfg.wp_timeout = 30;
cfg.r_sense = 1.2; cfg.safe_buffer = 0.06;  cfg.gamma = 2.0;
cfg.o0  = [1.3; 0.15];  cfg.vo0 = [-0.16; 0.16];  cfg.obs_r = 0.15;
cfg.obs_sa2_true = 0.0025;  cfg.obs_sa2_filter = 0.02;
cfg.track_drop = 3.0;

opt.N = 14; opt.T = 0.2; opt.rob_diam = 0.3;
opt.v_max = 0.6; opt.omega_max = pi/4;
opt.xy_min = -2; opt.xy_max = 2;
fprintf('Building parameterised-obstacle NMPC ...\n');
mpc = mc_build_mpc_dyn(opt);

L = size(cfg.lm,2);
maxiter = round(cfg.sim_tim/mpc.T);
ntot = maxiter + mpc.N + 2;

M = 50;

%% ---- EVALUATION (original paired seeds; frozen c from calibration) ----
L0 = load('calibration_result.mat'); c = L0.c;
fprintf('frozen c = %.4f m (from disjoint calibration stream)\n', c);
rng(2024);
noises = cell(M,1);
for kk = 1:M
    % original evaluation stream: sequential draws under rng(2024)
    nz.u = [sqrt(cfg.var_v); sqrt(cfg.var_w)] .* randn(2, maxiter);
    nz.z = zeros(2, L, maxiter);
    nz.z(1,:,:) = sqrt(cfg.var_d) * randn(1, L, maxiter);
    nz.z(2,:,:) = sqrt(cfg.var_a) * randn(1, L, maxiter);
    nz.oz   = [sqrt(cfg.var_d); sqrt(cfg.var_a)] .* randn(2, maxiter);
    nz.oacc = sqrt(cfg.obs_sa2_true) * randn(2, ntot);
    noises{kk} = nz;
end

arms = {'cv_cov', 0; 'cv_fixed', c};
R = struct();
for a = 1:2
    md = arms{a,1}; fe = arms{a,2};
    cc = cfg; if fe > 0, cc.fixed_extra = fe; end
    coll = false(M,1); mi_act = zeros(M,1); mi_all = zeros(M,1); afr = zeros(M,1);
    pk = zeros(M,1); md_i = zeros(M,1); iq = zeros(M,1);
    pathl = zeros(M,1); clr = zeros(M,1); sf = zeros(M,1); sm = zeros(M,1);
    for kk = 1:M
        r = mc_run_trial_um_dyn_v2(md, mpc, cc, noises{kk});
        coll(kk) = r.min_clear < 0; clr(kk) = r.min_clear;
        mi_act(kk) = r.mean_inflation_act; mi_all(kk) = r.mean_inflation_all;
        afr(kk) = r.active_frac; pk(kk) = r.peak_inflation;
        md_i(kk) = r.med_inflation_act; iq(kk) = r.iqr_inflation_act;
        pathl(kk) = r.path_len; sf(kk) = r.slack_freq; sm(kk) = r.slack_max;
    end
    R.(md) = struct('coll',coll,'clr',clr,'mi_act',mi_act,'mi_all',mi_all, ...
        'afr',afr,'pk',pk,'med',md_i,'iqr',iq,'path',pathl,'sf',sf,'sm',sm);
end
fid = fopen('dyn_ablation_v2_results.txt','w');
pr = @(varargin) fprintf(fid, varargin{:});
pr('=== Corrected dynamic matched-margin ablation (frozen protocol) ===\n');
pr('Comparator: fixed margin MATCHED TO THE REFERENCE ARM''S ACTIVE-TRACK MEAN\n');
pr('(not exposure-matched; residual mismatch reported below). c = %.4f m, frozen\n', c);
pr('from a disjoint calibration stream before any evaluation trial ran.\n\n');
for a = 1:2
    md = arms{a,1}; r = R.(md);
    nc = sum(r.coll);
    pr('%-9s coll %d/%d | active frac %.3f | infl|active mean %.3f med %.3f IQR %.3f max %.3f\n', ...
        md, nc, M, mean(r.afr), mean(r.mi_act), mean(r.med), mean(r.iqr), max(r.pk));
    pr('%-9s mission-integrated infl %.4f m | path %.2f | worst clear %+.3f | slack freq %.3f max %.3g\n', ...
        '', mean(r.mi_all), mean(r.path), min(r.clr), mean(r.sf), max(r.sm));
end
a1 = R.cv_cov.coll; a2 = R.cv_fixed.coll;
b10 = sum(a1 & ~a2); b01 = sum(~a1 & a2); nd = b10 + b01;
if nd > 0, pmc = min(1, 2*binocdf(min(b01,b10), nd, 0.5)); else, pmc = 1; end
pr('\npaired collisions: discordant cov-only %d : fixed-only %d | exact McNemar p = %.3g\n', b10, b01, pmc);
pr('risk difference %+.1f pp\n', 100*(sum(a2)-sum(a1))/M);
mm = 100*(mean(R.cv_fixed.mi_act) - mean(R.cv_cov.mi_act)) / mean(R.cv_cov.mi_act);
pr('\nresidual exposure mismatch (realised active-track means): %+.1f%%\n', mm);
if abs(mm) > 10
    pr('EXCEEDS the pre-declared 10%% bound: report as APPROXIMATELY matched;\n');
    pr('no causal size-vs-adaptivity conclusion may be drawn (frozen rule 5.3).\n');
end
pr('Pre-declared: a small discordant count is INCONCLUSIVE, not evidence of\n');
pr('either controller''s superiority (frozen rule 5.1).\n');
fclose(fid); type dyn_ablation_v2_results.txt
save('dyn_ablation_v2_results.mat','R','c','M');
fprintf('\nDYN_ABLATION_V2_DONE\n');
