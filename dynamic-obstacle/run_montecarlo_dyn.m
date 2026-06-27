% RUN_MONTECARLO_DYN  Dynamic-obstacle Monte-Carlo study. A robot (EKF-SLAM pose)
% must reach a goal while avoiding a MOVING obstacle tracked by a CV-EKF. Four
% obstacle-handling strategies are compared on identical noise realizations:
%   oracle   : true pose + clairvoyant true obstacle future (ceiling)
%   static   : freeze the obstacle at its current estimate (ignores motion)
%   cv_fixed : CV-predicted obstacle track, fixed safety margin
%   cv_cov   : CV-predicted track + covariance-aware margin (proposed)
%
% Robot localization (SLAM) is identical across static/cv_fixed/cv_cov, so the
% comparison isolates the obstacle-handling strategy.

clear; clc; close all;
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

% ---------------- scenario ----------------
cfg.lm = [-0.5 1 -1.5; 0.5 1 0.0];
cfg.xs = [1.5; 1.5; 0];  cfg.x0_nom = [0;0;0];
cfg.sigma_init_pos=0.10; cfg.sigma_init_th=0.05; cfg.lm_prior_std=0.10;
cfg.var_v=1e-6; cfg.var_w=1e-6; cfg.var_d=0.01; cfg.var_a=0.01;
cfg.sim_tim=30; cfg.tol=0.05; cfg.safe_buffer=0.06; cfg.gamma=2.0;
cfg.obs_r=0.15; cfg.obs_sa2_filter=0.02;     % tracker's assumed accel variance
cfg.o0=[1.3;0.15]; cfg.vo0=[-0.16;0.16];     % crossing obstacle
obs_sa_true=0.15;                             % true obstacle accel-noise std

opt.N=14; opt.T=0.2; opt.rob_diam=0.3; opt.v_max=0.6; opt.omega_max=pi/4; opt.rho=1e3;
fprintf('Building dynamic-obstacle NMPC ...\n');
mpc = mc_build_mpc_dyn(opt);
L=size(cfg.lm,2); maxiter=round(cfg.sim_tim/mpc.T); ntot=maxiter+opt.N+2;

M=50; modes={'oracle','static','cv_fixed','cv_cov'};
rng(2024);
noises=cell(M,1);
for k=1:M
    nz.offset=[cfg.sigma_init_pos*randn;cfg.sigma_init_pos*randn;cfg.sigma_init_th*randn];
    nz.u=[sqrt(cfg.var_v);sqrt(cfg.var_w)].*randn(2,maxiter);
    nz.z=zeros(2,L,maxiter); nz.z(1,:,:)=sqrt(cfg.var_d)*randn(1,L,maxiter); nz.z(2,:,:)=sqrt(cfg.var_a)*randn(1,L,maxiter);
    nz.lm=cfg.lm_prior_std*randn(2,L);
    nz.oacc=obs_sa_true*randn(2,ntot);
    nz.oz=zeros(2,maxiter); nz.oz(1,:)=sqrt(cfg.var_d)*randn(1,maxiter); nz.oz(2,:)=sqrt(cfg.var_a)*randn(1,maxiter);
    noises{k}=nz;
end

results=struct();
for mi=1:numel(modes)
    m=modes{mi}; COL=false(M,1); TR=zeros(M,1); MC=zeros(M,1); PL=zeros(M,1); sample=[];
    fprintf('Running %-9s ',m);
    for k=1:M
        r=mc_run_trial_dyn(m,mpc,cfg,noises{k});
        COL(k)=r.collided; TR(k)=r.term_true_ref; MC(k)=r.min_clear; PL(k)=r.path_len;
        if k==1, sample=r; end
    end
    results.(m).collided=COL; results.(m).term=TR; results.(m).minclear=MC;
    results.(m).pathlen=PL; results.(m).sample=sample;
    fprintf('done.\n');
end

ci=@(v) 1.96*std(v)/sqrt(numel(v));
fprintf('\n=== Dynamic-obstacle Monte-Carlo (M=%d, gamma=%.1f) ===\n',M,cfg.gamma);
fprintf('%-9s | %-11s | %-18s | %-14s | %-10s\n','mode','collision %','term true-ref [m]','min clear [m]','path [m]');
fprintf('%s\n',repmat('-',1,74));
for mi=1:numel(modes)
    m=modes{mi}; TR=results.(m).term;
    fprintf('%-9s | %8.1f    | %6.4f +/- %6.4f  | %+7.4f       | %6.3f\n', ...
        m, 100*mean(results.(m).collided), mean(TR), ci(TR), ...
        min(results.(m).minclear), mean(results.(m).pathlen));
end
fprintf('%s\n',repmat('-',1,74));
save('mc_results_dyn.mat','results','cfg','M');
fprintf('Saved mc_results_dyn.mat\n');

% ---- figures ----
% (1) collision rate
f1=figure('Color','w');
rates=cellfun(@(m) 100*mean(results.(m).collided), modes);
bar(rates,0.6); set(gca,'XTickLabel',modes,'FontName','Times New Roman','FontSize',12);
ylabel('Collision rate (%)'); title('Dynamic obstacle: collisions by strategy'); grid on
saveas(f1,'mc_dyn_collision.png');

% (2) example trajectories (trial 1): 4 robot paths + shared obstacle path
f2=figure('Color','w'); hold on
ob=results.oracle.sample.otraj;       % obstacle path is identical across modes (same noise)
plot(ob(1,:),ob(2,:),'--','Color',[0.6 0.2 0.2],'LineWidth',1.4,'DisplayName','obstacle path');
plot(ob(1,1),ob(2,1),'s','Color',[0.6 0.2 0.2],'MarkerFaceColor',[0.9 0.6 0.6]);
plot(ob(1,end),ob(2,end),'x','Color',[0.6 0.2 0.2],'MarkerSize',9);
plot(cfg.lm(1,:),cfg.lm(2,:),'k^','MarkerFaceColor','y','MarkerSize',8,'DisplayName','landmarks');
plot(cfg.xs(1),cfg.xs(2),'gp','MarkerFaceColor','g','MarkerSize',13,'DisplayName','goal');
cols={'k','m','c','b'};
for mi=1:numel(modes)
    tr=results.(modes{mi}).sample.traj;
    plot(tr(1,:),tr(2,:),cols{mi},'LineWidth',1.5,'DisplayName',modes{mi});
end
axis equal; grid on; legend('Location','eastoutside');
xlabel('x (m)'); ylabel('y (m)'); title('Dynamic obstacle: true trajectories (trial 1)');
set(gca,'FontName','Times New Roman','FontSize',12);
saveas(f2,'mc_dyn_trajectories.png');
fprintf('Saved mc_dyn_collision.png and mc_dyn_trajectories.png\n');
