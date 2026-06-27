% fig_side_by_side.m  Paper figure: static (collision) vs covariance-aware (safe)
% on the SAME trial (identical obstacle motion). Robot & obstacle drawn as disks
% at evenly spaced trajectory fractions; colour = progress; collision overlap
% highlighted. Saves PNG + vector PDF.
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

cfg.lm=[-0.5 1 -1.5;0.5 1 0]; cfg.xs=[1.5;1.5;0]; cfg.x0_nom=[0;0;0];
cfg.sigma_init_pos=0.10; cfg.sigma_init_th=0.05; cfg.lm_prior_std=0.10;
cfg.var_v=1e-6; cfg.var_w=1e-6; cfg.var_d=0.01; cfg.var_a=0.01;
cfg.sim_tim=30; cfg.tol=0.05; cfg.safe_buffer=0.06; cfg.gamma=2.0;
cfg.obs_r=0.15; cfg.obs_sa2_filter=0.02; cfg.o0=[1.3;0.15]; cfg.vo0=[-0.16;0.16];
obs_sa_true=0.15;
opt.N=14; opt.T=0.2; opt.rob_diam=0.3; opt.v_max=0.6; opt.omega_max=pi/4; opt.rho=1e3;
mpc=mc_build_mpc_dyn(opt); L=size(cfg.lm,2); maxiter=round(cfg.sim_tim/mpc.T); ntot=maxiter+opt.N+2;
rob_r=mpc.rob_r; sumr=rob_r+cfg.obs_r;

rng(2024);   % trial 1 = noises{1}
nz.offset=[cfg.sigma_init_pos*randn;cfg.sigma_init_pos*randn;cfg.sigma_init_th*randn];
nz.u=[sqrt(cfg.var_v);sqrt(cfg.var_w)].*randn(2,maxiter);
nz.z=zeros(2,L,maxiter); nz.z(1,:,:)=sqrt(cfg.var_d)*randn(1,L,maxiter); nz.z(2,:,:)=sqrt(cfg.var_a)*randn(1,L,maxiter);
nz.lm=cfg.lm_prior_std*randn(2,L);
nz.oacc=obs_sa_true*randn(2,ntot);
nz.oz=zeros(2,maxiter); nz.oz(1,:)=sqrt(cfg.var_d)*randn(1,maxiter); nz.oz(2,:)=sqrt(cfg.var_a)*randn(1,maxiter);

Rs = mc_run_trial_dyn('static', mpc, cfg, nz);
Rc = mc_run_trial_dyn('cv_cov', mpc, cfg, nz);
fprintf('static collided=%d (min %.3f) ; cv_cov collided=%d (min %.3f)\n', ...
        Rs.collided,Rs.min_clear,Rc.collided,Rc.min_clear);

% common axis limits over both trials
allx=[Rs.traj(1,:) Rs.otraj(1,:) Rc.traj(1,:) Rc.otraj(1,:) cfg.lm(1,:) cfg.xs(1)];
ally=[Rs.traj(2,:) Rs.otraj(2,:) Rc.traj(2,:) Rc.otraj(2,:) cfg.lm(2,:) cfg.xs(2)];
p=0.4; lims=[min(allx)-p max(allx)+p min(ally)-p max(ally)+p];

f=figure('Color','w','Position',[80 80 1180 560]);
t=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
draw_panel(nexttile, Rs, cfg, rob_r, sumr, lims, '(a) Static margin (ignores motion): collision');
draw_panel(nexttile, Rc, cfg, rob_r, sumr, lims, '(b) Covariance-aware margin: collision-free');
cb=colorbar; cb.Layout.Tile='east'; cb.Label.String='trajectory progress (%)';
colormap(parula); caxis([0 100]);

set(f,'PaperPositionMode','auto');
saveas(f,'fig_static_vs_cvcov.png');
print(f,'fig_static_vs_cvcov.pdf','-dpdf','-bestfit');
fprintf('Saved fig_static_vs_cvcov.png and .pdf\n');

% --------------------------------------------------------------
function draw_panel(ax, R, cfg, rob_r, sumr, lims, ttl)
hold(ax,'on');
tr=R.traj; ot=R.otraj; K=size(tr,2);
ang=linspace(0,2*pi,60); cxx=cos(ang); cyy=sin(ang);
clr=sqrt(sum((tr(1:2,:)-ot).^2,1))-sumr;
fracs=[0 0.2 0.4 0.6 0.8 1.0];
snaps=unique(max(1,round(fracs*(K-1))+1));
% faint paths
plot(ax,tr(1,:),tr(2,:),'-','Color',[0.55 0.55 1],'LineWidth',1.0);
plot(ax,ot(1,:),ot(2,:),'-','Color',[1 0.55 0.55],'LineWidth',1.0);
plot(ax,cfg.lm(1,:),cfg.lm(2,:),'k^','MarkerFaceColor','y','MarkerSize',8);
plot(ax,cfg.xs(1),cfg.xs(2),'gp','MarkerFaceColor','g','MarkerSize',15);
cmap=parula(256);
for kk=snaps
    fr=(kk-1)/max(K-1,1); col=cmap(max(1,round(fr*255)+1),:);
    fill(ax,tr(1,kk)+rob_r*cxx,tr(2,kk)+rob_r*cyy,col,'FaceAlpha',0.18,'EdgeColor',col,'LineWidth',1.4);
    th=tr(3,kk); plot(ax,tr(1,kk)+[0 rob_r*cos(th)],tr(2,kk)+[0 rob_r*sin(th)],'-','Color',col,'LineWidth',1.2);
    if clr(kk)<0, ofc=[0.95 0 0]; oew=2.4; else, ofc=[0.85 0.32 0.32]; oew=1.6; end
    fill(ax,ot(1,kk)+cfg.obs_r*cxx,ot(2,kk)+cfg.obs_r*cyy,ofc,'FaceAlpha',0.33,'EdgeColor',col,'LineWidth',oew);
end
% mark tightest approach
[mn,kc]=min(clr);
txt=sprintf('min clearance %+.2f m',mn);
text(ax,lims(1)+0.1,lims(4)-0.18,txt,'Color',(mn<0)*[0.7 0 0],'FontWeight','bold','FontSize',11);
axis(ax,'equal'); axis(ax,lims); grid(ax,'on'); box(ax,'on');
xlabel(ax,'x (m)'); ylabel(ax,'y (m)'); title(ax,ttl);
set(ax,'FontName','Times New Roman','FontSize',12);
end
