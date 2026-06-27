% draw_dyn_shapes.m  Visualize the moving obstacle and the robot as DISK SHAPES
% (not just paths) at several time instants of one trial.
load('mc_results_dyn.mat','results','cfg');
S  = results.cv_cov.sample;     % proposed method, trial 1
tr = S.traj;  ot = S.otraj;     % robot pose (3xK), obstacle pos (2xK)
K  = size(tr,2);
rob_r = 0.15;  obs_r = cfg.obs_r;
ang = linspace(0,2*pi,60);  cx = cos(ang);  cy = sin(ang);

nsnap = 7;
snaps = unique(round(linspace(1,K,nsnap)));
cmap  = parula(numel(snaps));

f = figure('Color','w','Position',[100 100 760 620]); hold on

% faint full paths
plot(tr(1,:),tr(2,:),'-','Color',[0.6 0.6 1],'LineWidth',1.0);
plot(ot(1,:),ot(2,:),'-','Color',[1 0.6 0.6],'LineWidth',1.0);

% landmarks + goal
plot(cfg.lm(1,:),cfg.lm(2,:),'k^','MarkerFaceColor','y','MarkerSize',9);
plot(cfg.xs(1),cfg.xs(2),'gp','MarkerFaceColor','g','MarkerSize',16);

for s = 1:numel(snaps)
    kk = snaps(s); col = cmap(s,:);
    % robot disk (filled, time-coloured outline)
    fill(tr(1,kk)+rob_r*cx, tr(2,kk)+rob_r*cy, col, 'FaceAlpha',0.20, ...
         'EdgeColor',col,'LineWidth',1.4);
    % robot heading tick
    th = tr(3,kk);
    plot(tr(1,kk)+[0 rob_r*cos(th)], tr(2,kk)+[0 rob_r*sin(th)],'-','Color',col,'LineWidth',1.4);
    % obstacle disk (filled red, same time colour outline)
    fill(ot(1,kk)+obs_r*cx, ot(2,kk)+obs_r*cy, [0.85 0.3 0.3], 'FaceAlpha',0.35, ...
         'EdgeColor',col,'LineWidth',1.6);
    % step labels
    text(tr(1,kk),tr(2,kk),sprintf('%d',kk),'FontSize',8,'HorizontalAlignment','center');
    text(ot(1,kk),ot(2,kk),sprintf('%d',kk),'FontSize',8,'HorizontalAlignment','center','Color','w');
end

% legend proxies
hr = fill(NaN,NaN,[0.4 0.4 1],'FaceAlpha',0.2,'EdgeColor',[0.2 0.2 1]);
ho = fill(NaN,NaN,[0.85 0.3 0.3],'FaceAlpha',0.35,'EdgeColor',[0.6 0.2 0.2]);
legend([hr ho],{'robot (r=0.15 m)','obstacle (r=0.15 m)'},'Location','northwest');

axis equal; grid on; box on
xlabel('x (m)'); ylabel('y (m)');
title('Robot and moving obstacle as disks (numbers = time step), cv\_cov trial 1');
cb = colorbar; cb.Label.String = 'time progression';
colormap(parula); caxis([snaps(1) snaps(end)]);
set(gca,'FontName','Times New Roman','FontSize',12);
saveas(f,'dyn_shapes_snapshots.png');
fprintf('Saved dyn_shapes_snapshots.png  (K=%d steps, snapshots at: %s)\n', K, mat2str(snaps));
