% RUN_GAMMA_SWEEP_DYN  Safety-vs-efficiency trade-off of the covariance-aware
% margin for the MOVING obstacle. Sweeps gamma for the cv_cov controller and
% records collision rate, detour (path length) and terminal accuracy. gamma=0
% reproduces the CV+fixed-margin baseline. static (ignore-motion) collision rate
% and the oracle path length are computed once as references.

clear; clc; close all;
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

cfg.lm=[-0.5 1 -1.5;0.5 1 0]; cfg.xs=[1.5;1.5;0]; cfg.x0_nom=[0;0;0];
cfg.sigma_init_pos=0.10; cfg.sigma_init_th=0.05; cfg.lm_prior_std=0.10;
cfg.var_v=1e-6; cfg.var_w=1e-6; cfg.var_d=0.01; cfg.var_a=0.01;
cfg.sim_tim=30; cfg.tol=0.05; cfg.safe_buffer=0.06; cfg.gamma=0;  % overridden in sweep
cfg.obs_r=0.15; cfg.obs_sa2_filter=0.02; cfg.o0=[1.3;0.15]; cfg.vo0=[-0.16;0.16];
obs_sa_true=0.15;
opt.N=14; opt.T=0.2; opt.rob_diam=0.3; opt.v_max=0.6; opt.omega_max=pi/4; opt.rho=1e3;
fprintf('Building NMPC ...\n');
mpc=mc_build_mpc_dyn(opt); L=size(cfg.lm,2); maxiter=round(cfg.sim_tim/mpc.T); ntot=maxiter+opt.N+2;

gammas=[0 0.5 1.0 1.5 2.0 2.5 3.0]; M=50;
rng(2024); noises=cell(M,1);
for k=1:M
    nz.offset=[cfg.sigma_init_pos*randn;cfg.sigma_init_pos*randn;cfg.sigma_init_th*randn];
    nz.u=[sqrt(cfg.var_v);sqrt(cfg.var_w)].*randn(2,maxiter);
    nz.z=zeros(2,L,maxiter); nz.z(1,:,:)=sqrt(cfg.var_d)*randn(1,L,maxiter); nz.z(2,:,:)=sqrt(cfg.var_a)*randn(1,L,maxiter);
    nz.lm=cfg.lm_prior_std*randn(2,L);
    nz.oacc=obs_sa_true*randn(2,ntot);
    nz.oz=zeros(2,maxiter); nz.oz(1,:)=sqrt(cfg.var_d)*randn(1,maxiter); nz.oz(2,:)=sqrt(cfg.var_a)*randn(1,maxiter);
    noises{k}=nz;
end
plen=@(tr) sum(sqrt(sum(diff(tr(1:2,:),1,2).^2,1)));

% references
colS=false(M,1); lenO=zeros(M,1);
fprintf('References (static, oracle) ...\n');
for k=1:M
    rs=mc_run_trial_dyn('static',mpc,cfg,noises{k}); colS(k)=rs.collided;
    ro=mc_run_trial_dyn('oracle',mpc,cfg,noises{k}); lenO(k)=plen(ro.traj);
end
static_rate=100*mean(colS); oracle_len=mean(lenO);

G=numel(gammas); col=zeros(G,1); pth=zeros(G,1); ter=zeros(G,1);
for gi=1:G
    c=cfg; c.gamma=gammas(gi);
    COL=false(M,1); PL=zeros(M,1); TR=zeros(M,1);
    fprintf('gamma=%.1f ... ',gammas(gi));
    for k=1:M
        r=mc_run_trial_dyn('cv_cov',mpc,c,noises{k});
        COL(k)=r.collided; PL(k)=plen(r.traj); TR(k)=r.term_true_ref;
    end
    col(gi)=100*mean(COL); pth(gi)=mean(PL); ter(gi)=mean(TR);
    fprintf('collision=%4.1f%%  path=%.3f m  term=%.3f m\n',col(gi),pth(gi),ter(gi));
end

fprintf('\n=== Dynamic gamma sweep (M=%d). Ref: static collision %.1f%%, oracle path %.3f m ===\n',M,static_rate,oracle_len);
fprintf('%-6s | %-12s | %-12s | %-12s\n','gamma','collision %','path [m]','term err [m]');
fprintf('%s\n',repmat('-',1,52));
for gi=1:G, fprintf('%-6.1f | %10.1f   | %10.3f   | %10.3f\n',gammas(gi),col(gi),pth(gi),ter(gi)); end
fprintf('%s\n',repmat('-',1,52));
save('gamma_sweep_dyn.mat','gammas','col','pth','ter','static_rate','oracle_len','M');

f=figure('Color','w');
yyaxis left
plot(gammas,col,'-o','LineWidth',1.8); hold on
yline(static_rate,'--','static (ignore motion)','Color',[0.7 0 0],'LineWidth',1.2,'LabelHorizontalAlignment','left');
ylabel('Collision rate (%)'); ylim([-2 max(static_rate,max(col))+4]);
yyaxis right
plot(gammas,pth,'-s','LineWidth',1.8);
yline(oracle_len,':','oracle path','LineWidth',1.2,'LabelHorizontalAlignment','right');
ylabel('Mean path length (m)');
xlabel('\gamma (chance-constraint factor)');
title('Dynamic obstacle: safety vs. efficiency trade-off');
grid on; set(gca,'FontName','Times New Roman','FontSize',12);
legend({'cv\_cov collision rate','static (ref)','cv\_cov path length','oracle path'},'Location','east');
saveas(f,'gamma_sweep_dyn.png');
fprintf('Saved gamma_sweep_dyn.mat and gamma_sweep_dyn.png\n');
