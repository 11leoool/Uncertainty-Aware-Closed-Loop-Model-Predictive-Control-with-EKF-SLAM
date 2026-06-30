% run_mc100.m  Re-run ALL Monte-Carlo studies at M=100 to tighten the collision
% confidence intervals (Wilson 95%). Prints a consolidated summary.
% (All local functions are defined at the END, as MATLAB requires.)
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*
M = 100;
fprintf('==================  MONTE CARLO @ M=%d  ==================\n', M);

base.lm=[-0.5 1 -1.5;0.5 1 0]; base.xs=[1.5;1.5;0]; base.x0_nom=[0;0;0];
base.sigma_init_pos=0.10; base.sigma_init_th=0.05; base.lm_prior_std=0.10;
base.var_v=1e-6; base.var_w=1e-6; base.var_d=0.01; base.var_a=0.01; base.tol=0.05;
L=size(base.lm,2);
ci = @(v) 1.96*std(v)/sqrt(numel(v));

% =================== 1. FREE-SPACE TRACKING ===================
cfg=base; cfg.sim_tim=20; mpc=mc_build_mpc(); mi=round(cfg.sim_tim/mpc.T);
rng(2024); no=cell(M,1); for k=1:M, no{k}=noise_obs(cfg,L,mi); end
modes={'oracle','odom','slam'}; fprintf('\n[1] Free-space tracking (M=%d)\n',M);
for j=1:3, TR=zeros(M,1); LE=zeros(M,1);
  for k=1:M, r=mc_run_trial(modes{j},mpc,cfg,no{k}); TR(k)=r.term_true_ref; LE(k)=r.mean_true_est; end
  fprintf('  %-7s term=%.4f +/- %.4f  loc=%.4f +/- %.4f\n',modes{j},mean(TR),ci(TR),mean(LE),ci(LE)); end

% =================== obstacle setup ===================
opt.N=14; opt.T=0.2; opt.rob_diam=0.3; opt.v_max=0.6; opt.omega_max=pi/4;
opt.obs=[0.5 0.5 0.15]; opt.xy_min=-2; opt.xy_max=2;
mpcO=mc_build_mpc_obs(opt); cfgO=base; cfgO.sim_tim=30; cfgO.safe_buffer=0.06; cfgO.gamma=2.0;
miO=round(cfgO.sim_tim/opt.T);
rng(2024); noO=cell(M,1); for k=1:M, noO{k}=noise_obs(cfgO,L,miO); end

fprintf('\n[2] Static obstacle, Stage A (fixed margin), M=%d\n',M);
for j=1:3, COL=false(M,1); TR=zeros(M,1); c=cfgO; c.cov_aware=false;
  for k=1:M, r=mc_run_trial_obs(modes{j},mpcO,c,noO{k}); COL(k)=r.collided; TR(k)=r.term_true_ref; end
  fprintf('  %-7s collision %s  term=%.4f\n',modes{j},wil(sum(COL),M),mean(TR)); end
fprintf('[3] Static obstacle, Stage B (covariance-aware), M=%d\n',M);
c=cfgO; c.cov_aware=true; COL=false(M,1); TR=zeros(M,1);
for k=1:M, r=mc_run_trial_obs('slam',mpcO,c,noO{k}); COL(k)=r.collided; TR(k)=r.term_true_ref; end
fprintf('  slam(cov) collision %s  term=%.4f\n',wil(sum(COL),M),mean(TR));

fprintf('[4] Static gamma sweep (slam cov-aware), M=%d\n',M);
for g=[0 0.5 1 1.5 2 2.5 3], c=cfgO; c.cov_aware=true; c.gamma=g; COL=false(M,1); PL=zeros(M,1);
  for k=1:M, r=mc_run_trial_obs('slam',mpcO,c,noO{k}); COL(k)=r.collided; PL(k)=r.path_len; end
  fprintf('  gamma=%.1f collision %s  path=%.3f\n',g,wil(sum(COL),M),mean(PL)); end

% =================== dynamic setup ===================
optd=opt; optd.rho=1e3; mpcD=mc_build_mpc_dyn(optd);
cfgD=cfgO; cfgD.obs_r=0.15; cfgD.obs_sa2_filter=0.02; cfgD.o0=[1.3;0.15]; cfgD.vo0=[-0.16;0.16];
sat=0.15; miD=round(cfgD.sim_tim/optd.T); ntot=miD+optd.N+2;
rng(2024); noD=cell(M,1); for k=1:M, noD{k}=noise_dyn(cfgD,L,miD,ntot,sat); end

fprintf('\n[5] Dynamic obstacle, 4 strategies, M=%d\n',M);
dm={'oracle','static','cv_fixed','cv_cov'};
for j=1:4, COL=false(M,1); TR=zeros(M,1); PL=zeros(M,1); c=cfgD; c.gamma=2.0;
  for k=1:M, r=mc_run_trial_dyn(dm{j},mpcD,c,noD{k}); COL(k)=r.collided; TR(k)=r.term_true_ref; PL(k)=r.path_len; end
  fprintf('  %-8s collision %s  term=%.4f  path=%.3f\n',dm{j},wil(sum(COL),M),mean(TR),mean(PL)); end

fprintf('[6] Dynamic gamma sweep (cv_cov), M=%d\n',M);
for g=[0 0.5 1 1.5 2 2.5 3], c=cfgD; c.gamma=g; COL=false(M,1); PL=zeros(M,1);
  for k=1:M, r=mc_run_trial_dyn('cv_cov',mpcD,c,noD{k}); COL(k)=r.collided; PL(k)=r.path_len; end
  fprintf('  gamma=%.1f collision %s  path=%.3f\n',g,wil(sum(COL),M),mean(PL)); end

% =================== 7. ABLATION ===================
fprintf('\n[7] Ablation: margin size vs adaptivity, M=%d\n',M);
c=cfgO; c.cov_aware=true; c.gamma=2; col1=false(M,1); p1=zeros(M,1); inf1=zeros(M,1); st1=zeros(M,1);
for k=1:M, r=mc_run_trial_obs('slam',mpcO,c,noO{k}); col1(k)=r.collided; p1(k)=r.path_len; inf1(k)=r.infl_mean; st1(k)=r.steps; end
cmS=sum(inf1.*st1)/sum(st1);
cF=cfgO; cF.cov_aware=false; cF.safe_buffer=0.06+cmS; col2=false(M,1); p2=zeros(M,1);
for k=1:M, r=mc_run_trial_obs('slam',mpcO,cF,noO{k}); col2(k)=r.collided; p2(k)=r.path_len; end
cB=cfgO; cB.cov_aware=false; col0=false(M,1); p0=zeros(M,1);
for k=1:M, r=mc_run_trial_obs('slam',mpcO,cB,noO{k}); col0(k)=r.collided; p0(k)=r.path_len; end
fprintf('  STATIC  (matched extra=%.3f m)\n',cmS);
fprintf('    cv_fixed       collision %s  path=%.3f\n',wil(sum(col0),M),mean(p0));
fprintf('    fixed-matched  collision %s  path=%.3f\n',wil(sum(col2),M),mean(p2));
fprintf('    cv_cov         collision %s  path=%.3f\n',wil(sum(col1),M),mean(p1));
c=cfgD; c.gamma=2; dcol1=false(M,1); dp1=zeros(M,1); dinf=zeros(M,1); dst=zeros(M,1);
for k=1:M, r=mc_run_trial_dyn('cv_cov',mpcD,c,noD{k}); dcol1(k)=r.collided; dp1(k)=r.path_len; dinf(k)=r.infl_mean; dst(k)=r.steps; end
cmD=sum(dinf.*dst)/sum(dst);
cF=cfgD; cF.fixed_extra=cmD; dcol2=false(M,1); dp2=zeros(M,1);
for k=1:M, r=mc_run_trial_dyn('cv_fixedmatch',mpcD,cF,noD{k}); dcol2(k)=r.collided; dp2(k)=r.path_len; end
dcol0=false(M,1); dp0=zeros(M,1);
for k=1:M, r=mc_run_trial_dyn('cv_fixed',mpcD,cfgD,noD{k}); dcol0(k)=r.collided; dp0(k)=r.path_len; end
fprintf('  DYNAMIC (matched extra=%.3f m)\n',cmD);
fprintf('    cv_fixed       collision %s  path=%.3f\n',wil(sum(dcol0),M),mean(dp0));
fprintf('    fixed-matched  collision %s  path=%.3f\n',wil(sum(dcol2),M),mean(dp2));
fprintf('    cv_cov         collision %s  path=%.3f\n',wil(sum(dcol1),M),mean(dp1));
fprintf('\n==================  MC@%d DONE  ==================\n',M);

% =================== local functions (must be at end) ===================
function s=wil(x,n)
  p=x/n; z=1.96; d=1+z^2/n; c=(p+z^2/(2*n))/d;
  h=z*sqrt(p*(1-p)/n+z^2/(4*n^2))/d;
  s=sprintf('%4.1f%% [%4.1f, %4.1f]',100*p,100*max(0,c-h),100*min(1,c+h));
end
function nz=noise_obs(cfg,L,mi)
  nz.offset=[cfg.sigma_init_pos*randn;cfg.sigma_init_pos*randn;cfg.sigma_init_th*randn];
  nz.u=[sqrt(cfg.var_v);sqrt(cfg.var_w)].*randn(2,mi);
  nz.z=zeros(2,L,mi); nz.z(1,:,:)=sqrt(cfg.var_d)*randn(1,L,mi); nz.z(2,:,:)=sqrt(cfg.var_a)*randn(1,L,mi);
  nz.lm=cfg.lm_prior_std*randn(2,L);
end
function nz=noise_dyn(cfg,L,mi,ntot,sat)
  nz=noise_obs(cfg,L,mi);
  nz.oacc=sat*randn(2,ntot);
  nz.oz=zeros(2,mi); nz.oz(1,:)=sqrt(cfg.var_d)*randn(1,mi); nz.oz(2,:)=sqrt(cfg.var_a)*randn(1,mi);
end
