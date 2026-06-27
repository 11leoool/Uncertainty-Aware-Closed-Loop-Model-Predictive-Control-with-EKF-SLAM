% make_dyn_media.m  Build animations (GIF) and a collision snapshot showing the
% robot and moving obstacle as DISKS. Re-runs trials with the same seeds as the
% driver, captures a safe cv_cov run and the first colliding 'static' run.
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

% ---- config identical to run_montecarlo_dyn.m ----
cfg.lm=[-0.5 1 -1.5;0.5 1 0]; cfg.xs=[1.5;1.5;0]; cfg.x0_nom=[0;0;0];
cfg.sigma_init_pos=0.10; cfg.sigma_init_th=0.05; cfg.lm_prior_std=0.10;
cfg.var_v=1e-6; cfg.var_w=1e-6; cfg.var_d=0.01; cfg.var_a=0.01;
cfg.sim_tim=30; cfg.tol=0.05; cfg.safe_buffer=0.06; cfg.gamma=2.0;
cfg.obs_r=0.15; cfg.obs_sa2_filter=0.02; cfg.o0=[1.3;0.15]; cfg.vo0=[-0.16;0.16];
obs_sa_true=0.15;
opt.N=14; opt.T=0.2; opt.rob_diam=0.3; opt.v_max=0.6; opt.omega_max=pi/4; opt.rho=1e3;
mpc=mc_build_mpc_dyn(opt); L=size(cfg.lm,2); maxiter=round(cfg.sim_tim/mpc.T); ntot=maxiter+opt.N+2;
rob_r=mpc.rob_r; sumr=rob_r+cfg.obs_r;

M=50; rng(2024); noises=cell(M,1);
for k=1:M
    nz.offset=[cfg.sigma_init_pos*randn;cfg.sigma_init_pos*randn;cfg.sigma_init_th*randn];
    nz.u=[sqrt(cfg.var_v);sqrt(cfg.var_w)].*randn(2,maxiter);
    nz.z=zeros(2,L,maxiter); nz.z(1,:,:)=sqrt(cfg.var_d)*randn(1,L,maxiter); nz.z(2,:,:)=sqrt(cfg.var_a)*randn(1,L,maxiter);
    nz.lm=cfg.lm_prior_std*randn(2,L);
    nz.oacc=obs_sa_true*randn(2,ntot);
    nz.oz=zeros(2,maxiter); nz.oz(1,:)=sqrt(cfg.var_d)*randn(1,maxiter); nz.oz(2,:)=sqrt(cfg.var_a)*randn(1,maxiter);
    noises{k}=nz;
end

Rsafe = mc_run_trial_dyn('cv_cov', mpc, cfg, noises{1});
Rcol = []; colk = 0;
for k=1:M
    R = mc_run_trial_dyn('static', mpc, cfg, noises{k});
    if R.collided, Rcol=R; colk=k; break; end
end
fprintf('cv_cov trial1 collided=%d ; first static collision at trial %d\n', Rsafe.collided, colk);

make_gif(Rsafe, cfg, rob_r, sumr, 'dyn_anim_cvcov_safe.gif', 'cv\_cov (safe)');
if ~isempty(Rcol)
    make_gif(Rcol, cfg, rob_r, sumr, 'dyn_anim_static_collision.gif', 'static (collision)');
    snapshot_png(Rcol, cfg, rob_r, sumr, 'dyn_collision_snapshots.png');
end
fprintf('MEDIA_DONE\n');

% ============================ helpers ============================
function lims = data_lims(tr, ot, cfg)
ax=[tr(1,:) ot(1,:) cfg.lm(1,:) cfg.xs(1)]; ay=[tr(2,:) ot(2,:) cfg.lm(2,:) cfg.xs(2)];
p=0.45; lims=[min(ax)-p max(ax)+p min(ay)-p max(ay)+p];
end

function make_gif(R, cfg, rob_r, sumr, fname, ttl)
tr=R.traj; ot=R.otraj; K=size(tr,2); ang=linspace(0,2*pi,50); cxx=cos(ang); cyy=sin(ang);
lims=data_lims(tr,ot,cfg);
fig=figure('visible','off','Color','w','Position',[100 100 640 560]);
for k=1:K
    clf; hold on
    plot(cfg.lm(1,:),cfg.lm(2,:),'k^','MarkerFaceColor','y','MarkerSize',8);
    plot(cfg.xs(1),cfg.xs(2),'gp','MarkerFaceColor','g','MarkerSize',15);
    plot(tr(1,1:k),tr(2,1:k),'-','Color',[0.5 0.5 1],'LineWidth',1.0);
    plot(ot(1,1:k),ot(2,1:k),'-','Color',[1 0.5 0.5],'LineWidth',1.0);
    clr=norm(tr(1:2,k)-ot(:,k))-sumr;
    if clr<0, oc=[0.9 0 0]; else, oc=[0.85 0.35 0.35]; end
    fill(tr(1,k)+rob_r*cxx, tr(2,k)+rob_r*cyy,[0.3 0.5 1],'FaceAlpha',0.55,'EdgeColor',[0 0 0.6],'LineWidth',1.2);
    th=tr(3,k); plot(tr(1,k)+[0 rob_r*cos(th)],tr(2,k)+[0 rob_r*sin(th)],'k-','LineWidth',1.3);
    fill(ot(1,k)+cfg.obs_r*cxx, ot(2,k)+cfg.obs_r*cyy, oc,'FaceAlpha',0.6,'EdgeColor',[0.5 0 0],'LineWidth',1.2);
    axis equal; axis(lims); grid on; box on
    xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('%s   step %d/%d   clearance = %+.2f m', ttl, k, K, clr));
    set(gca,'FontName','Times New Roman','FontSize',11);
    drawnow; fr=getframe(fig); im=frame2im(fr); [A,map]=rgb2ind(im,256);
    if k==1, imwrite(A,map,fname,'gif','LoopCount',Inf,'DelayTime',0.12);
    else,    imwrite(A,map,fname,'gif','WriteMode','append','DelayTime',0.12); end
end
close(fig); fprintf('Saved %s (%d frames)\n', fname, K);
end

function snapshot_png(R, cfg, rob_r, sumr, fname)
tr=R.traj; ot=R.otraj; K=size(tr,2); ang=linspace(0,2*pi,60); cxx=cos(ang); cyy=sin(ang);
clr=sqrt(sum((tr(1:2,:)-ot).^2,1))-sumr; [~,kc]=min(clr);   % tightest (collision) step
snaps=unique([round(linspace(1,K,6)) kc]); cmap=parula(numel(snaps));
f=figure('Color','w','Position',[100 100 760 620]); hold on
plot(tr(1,:),tr(2,:),'-','Color',[0.6 0.6 1],'LineWidth',1.0);
plot(ot(1,:),ot(2,:),'-','Color',[1 0.6 0.6],'LineWidth',1.0);
plot(cfg.lm(1,:),cfg.lm(2,:),'k^','MarkerFaceColor','y','MarkerSize',9);
plot(cfg.xs(1),cfg.xs(2),'gp','MarkerFaceColor','g','MarkerSize',16);
for s=1:numel(snaps)
    kk=snaps(s); col=cmap(s,:); iscol = (clr(kk)<0);
    fill(tr(1,kk)+rob_r*cxx,tr(2,kk)+rob_r*cyy,col,'FaceAlpha',0.20,'EdgeColor',col,'LineWidth',1.4);
    th=tr(3,kk); plot(tr(1,kk)+[0 rob_r*cos(th)],tr(2,kk)+[0 rob_r*sin(th)],'-','Color',col,'LineWidth',1.3);
    if iscol, ofc=[0.95 0 0]; oew=2.2; else, ofc=[0.85 0.3 0.3]; oew=1.6; end
    fill(ot(1,kk)+cfg.obs_r*cxx,ot(2,kk)+cfg.obs_r*cyy,ofc,'FaceAlpha',0.35,'EdgeColor',col,'LineWidth',oew);
    text(tr(1,kk),tr(2,kk),sprintf('%d',kk),'FontSize',8,'HorizontalAlignment','center');
    text(ot(1,kk),ot(2,kk),sprintf('%d',kk),'FontSize',8,'HorizontalAlignment','center','Color','w');
end
text(ot(1,kc),ot(2,kc)-0.28,sprintf('collision (step %d, %.2f m)',kc,clr(kc)), ...
     'Color',[0.7 0 0],'FontWeight','bold','HorizontalAlignment','center');
axis equal; grid on; box on; xlabel('x (m)'); ylabel('y (m)');
title('static baseline: robot and obstacle disks overlap at collision');
cb=colorbar; cb.Label.String='time progression'; colormap(parula); caxis([snaps(1) snaps(end)]);
set(gca,'FontName','Times New Roman','FontSize',12);
saveas(f,fname); fprintf('Saved %s (collision step %d)\n', fname, kc);
end
