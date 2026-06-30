% RUN_ABLATION  Adaptivity-vs-size ablation (Devil's Advocate experiment).
% Question: is the covariance-aware safety gain due to ADAPTING the margin to live
% uncertainty, or merely to a LARGER average margin? Control = a FIXED margin set to
% cv_cov's MEAN applied inflation. Same mean margin, different shaping:
%   - if fixed-matched is as safe at equal/shorter path -> "size" explains it
%   - if cv_cov is safer per unit detour            -> "adaptivity" explains it
% Run for STATIC and DYNAMIC obstacles, same seeds as the main studies.

clear; clc; close all;
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

wil = @(x,n) deal_wilson(x,n);   % placeholder (defined below as local fn)

% ============================ STATIC ============================
cfg.lm=[-0.5 1 -1.5;0.5 1 0]; cfg.xs=[1.5;1.5;0]; cfg.x0_nom=[0;0;0];
cfg.sigma_init_pos=0.10; cfg.sigma_init_th=0.05; cfg.lm_prior_std=0.10;
cfg.var_v=1e-6; cfg.var_w=1e-6; cfg.var_d=0.01; cfg.var_a=0.01;
cfg.sim_tim=30; cfg.tol=0.05; cfg.safe_buffer=0.06; cfg.gamma=2.0; cfg.cov_aware=false;
opt.N=14; opt.T=0.2; opt.rob_diam=0.3; opt.v_max=0.6; opt.omega_max=pi/4;
opt.obs=[0.5 0.5 0.15]; opt.xy_min=-2; opt.xy_max=2;
fprintf('Building static NMPC ...\n'); mpcS=mc_build_mpc_obs(opt);
L=size(cfg.lm,2); maxiterS=round(cfg.sim_tim/opt.T);

M=50; rng(2024); noiseS=cell(M,1);
for k=1:M
  nz.offset=[cfg.sigma_init_pos*randn;cfg.sigma_init_pos*randn;cfg.sigma_init_th*randn];
  nz.u=[sqrt(cfg.var_v);sqrt(cfg.var_w)].*randn(2,maxiterS);
  nz.z=zeros(2,L,maxiterS); nz.z(1,:,:)=sqrt(cfg.var_d)*randn(1,L,maxiterS); nz.z(2,:,:)=sqrt(cfg.var_a)*randn(1,L,maxiterS);
  nz.lm=cfg.lm_prior_std*randn(2,L); noiseS{k}=nz;
end
% pass 1: cv_cov -> collisions, path, mean inflation
cS=cfg; cS.cov_aware=true; cS.gamma=2.0;
[colCov,pthCov,inflCov,stp]=deal(zeros(M,1),zeros(M,1),zeros(M,1),zeros(M,1));
for k=1:M, r=mc_run_trial_obs('slam',mpcS,cS,noiseS{k}); colCov(k)=r.collided; pthCov(k)=r.path_len; inflCov(k)=r.infl_mean; stp(k)=r.steps; end
cmatchS = sum(inflCov.*stp)/sum(stp);          % steps-weighted mean inflation
% pass 2: fixed margin matched to that mean
cF=cfg; cF.cov_aware=false; cF.safe_buffer=0.06+cmatchS;
[colFix,pthFix]=deal(zeros(M,1),zeros(M,1));
for k=1:M, r=mc_run_trial_obs('slam',mpcS,cF,noiseS{k}); colFix(k)=r.collided; pthFix(k)=r.path_len; end
% reference: plain fixed buffer (gamma=0)
cB=cfg; cB.cov_aware=false;
[colB,pthB]=deal(zeros(M,1),zeros(M,1));
for k=1:M, r=mc_run_trial_obs('slam',mpcS,cB,noiseS{k}); colB(k)=r.collided; pthB(k)=r.path_len; end

% ============================ DYNAMIC ============================
dcfg=cfg; dcfg.obs_r=0.15; dcfg.obs_sa2_filter=0.02; dcfg.o0=[1.3;0.15]; dcfg.vo0=[-0.16;0.16];
obs_sa_true=0.15;
optd=opt; optd.rho=1e3;
fprintf('Building dynamic NMPC ...\n'); mpcD=mc_build_mpc_dyn(optd);
maxiterD=round(dcfg.sim_tim/optd.T); ntot=maxiterD+optd.N+2;
rng(2024); noiseD=cell(M,1);
for k=1:M
  nz.offset=[dcfg.sigma_init_pos*randn;dcfg.sigma_init_pos*randn;dcfg.sigma_init_th*randn];
  nz.u=[sqrt(dcfg.var_v);sqrt(dcfg.var_w)].*randn(2,maxiterD);
  nz.z=zeros(2,L,maxiterD); nz.z(1,:,:)=sqrt(dcfg.var_d)*randn(1,L,maxiterD); nz.z(2,:,:)=sqrt(dcfg.var_a)*randn(1,L,maxiterD);
  nz.lm=dcfg.lm_prior_std*randn(2,L);
  nz.oacc=obs_sa_true*randn(2,ntot);
  nz.oz=zeros(2,maxiterD); nz.oz(1,:)=sqrt(dcfg.var_d)*randn(1,maxiterD); nz.oz(2,:)=sqrt(dcfg.var_a)*randn(1,maxiterD);
  noiseD{k}=nz;
end
% pass 1: cv_cov
dC=dcfg; dC.gamma=2.0;
[dcolCov,dpthCov,dinfl,dstp]=deal(zeros(M,1),zeros(M,1),zeros(M,1),zeros(M,1));
for k=1:M, r=mc_run_trial_dyn('cv_cov',mpcD,dC,noiseD{k}); dcolCov(k)=r.collided; dpthCov(k)=r.path_len; dinfl(k)=r.infl_mean; dstp(k)=r.steps; end
cmatchD = sum(dinfl.*dstp)/sum(dstp);
% pass 2: cv prediction + fixed margin matched to mean
dF=dcfg; dF.fixed_extra=cmatchD;
[dcolFix,dpthFix]=deal(zeros(M,1),zeros(M,1));
for k=1:M, r=mc_run_trial_dyn('cv_fixedmatch',mpcD,dF,noiseD{k}); dcolFix(k)=r.collided; dpthFix(k)=r.path_len; end
% reference: cv_fixed (delta0 only)
[dcolB,dpthB]=deal(zeros(M,1),zeros(M,1));
for k=1:M, r=mc_run_trial_dyn('cv_fixed',mpcD,dcfg,noiseD{k}); dcolB(k)=r.collided; dpthB(k)=r.path_len; end

% ============================ REPORT ============================
prn=@(name,col,pth) report_row(name,col,pth);
fprintf('\n================ ADAPTIVITY-VS-SIZE ABLATION (M=%d) ================\n',M);
fprintf('\n--- STATIC obstacle ---  cv_cov mean inflation matched as fixed extra = %.4f m\n',cmatchS);
fprintf('%-22s | %-22s | %-14s\n','controller','collision % [95%% Wilson]','mean path [m]');
fprintf('%s\n',repmat('-',1,64));
prn('cv_fixed (delta0)',colB,pthB);
prn('fixed-matched (size)',colFix,pthFix);
prn('cv_cov (adaptive)',colCov,pthCov);

fprintf('\n--- DYNAMIC obstacle ---  cv_cov mean inflation matched as fixed extra = %.4f m\n',cmatchD);
fprintf('%-22s | %-22s | %-14s\n','controller','collision % [95%% Wilson]','mean path [m]');
fprintf('%s\n',repmat('-',1,64));
prn('cv_fixed (delta0)',dcolB,dpthB);
prn('fixed-matched (size)',dcolFix,dpthFix);
prn('cv_cov (adaptive)',dcolCov,dpthCov);

save('ablation_results.mat','colB','colFix','colCov','pthB','pthFix','pthCov','cmatchS',...
     'dcolB','dcolFix','dcolCov','dpthB','dpthFix','dpthCov','cmatchD','M');
fprintf('\nSaved ablation_results.mat\n');

% bar figure: collision% and path, cv_cov vs fixed-matched, static + dynamic
fig=figure('Color','w','Position',[100 100 820 340]);
subplot(1,2,1);
bar([100*mean(colFix) 100*mean(colCov); 100*mean(dcolFix) 100*mean(dcolCov)]);
set(gca,'XTickLabel',{'static','dynamic'}); ylabel('Collision rate (%)');
legend({'fixed-matched','cv\_cov (adaptive)'},'Location','best'); title('Same mean margin: collisions'); grid on
subplot(1,2,2);
bar([mean(pthFix) mean(pthCov); mean(dpthFix) mean(dpthCov)]);
set(gca,'XTickLabel',{'static','dynamic'}); ylabel('Mean path length (m)');
legend({'fixed-matched','cv\_cov (adaptive)'},'Location','best'); title('Same mean margin: path'); grid on
saveas(fig,'ablation_compare.png');
fprintf('Saved ablation_compare.png\n');

% ---------------- local functions ----------------
function report_row(name,col,pth)
  x=sum(col); n=numel(col); p=x/n; z=1.96; d=1+z^2/n;
  c=(p+z^2/(2*n))/d; h=z*sqrt(p*(1-p)/n+z^2/(4*n^2))/d;
  fprintf('%-22s | %5.1f  [%4.1f, %4.1f]      | %6.3f\n',name,100*p,100*max(0,c-h),100*min(1,c+h),mean(pth));
end
function varargout=deal_wilson(varargin), varargout{1}=[]; end
