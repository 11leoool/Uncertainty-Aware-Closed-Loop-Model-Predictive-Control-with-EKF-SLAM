% time_perf_obs.m  Per-step timing of the obstacle-avoidance NMPC + EKF-SLAM loop
% (multiple shooting, N=14, T=0.2, one disk obstacle). slam mode.
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5'); import casadi.*

cfg.lm=[-0.5 1 -1.5;0.5 1 0]; cfg.xs=[1.5;1.5;0]; cfg.x0_nom=[0;0;0];
cfg.sigma_init_pos=0.10; cfg.sigma_init_th=0.05; cfg.lm_prior_std=0.10;
cfg.var_v=1e-6; cfg.var_w=1e-6; cfg.var_d=0.01; cfg.var_a=0.01;
cfg.sim_tim=30; cfg.tol=0.05; cfg.safe_buffer=0.06; cfg.cov_aware=false;
opt.N=14; opt.T=0.2; opt.rob_diam=0.3; opt.v_max=0.6; opt.omega_max=pi/4;
opt.obs=[0.5 0.5 0.15]; opt.xy_min=-2; opt.xy_max=2;
mpc=mc_build_mpc_obs(opt); L=size(cfg.lm,2); maxiter=round(cfg.sim_tim/mpc.T);
ns=mpc.n_states; nc=mpc.n_controls; N=mpc.N; T=mpc.T;

rng(2024); t_solve=[]; t_ekf=[]; nTrials=20;
for tr=1:nTrials
  nz.offset=[cfg.sigma_init_pos*randn;cfg.sigma_init_pos*randn;cfg.sigma_init_th*randn];
  nz.u=[sqrt(cfg.var_v);sqrt(cfg.var_w)].*randn(2,maxiter);
  nz.z=zeros(2,L,maxiter); nz.z(1,:,:)=sqrt(cfg.var_d)*randn(1,L,maxiter); nz.z(2,:,:)=sqrt(cfg.var_a)*randn(1,L,maxiter);
  nz.lm=cfg.lm_prior_std*randn(2,L);
  n=3+2*L; X=zeros(n,1); X(1:3)=cfg.x0_nom; X(4:end)=cfg.lm(:)+nz.lm(:);
  Sigma=zeros(n); Sigma(1:3,1:3)=diag([cfg.sigma_init_pos^2,cfg.sigma_init_pos^2,cfg.sigma_init_th^2]); Sigma(4:end,4:end)=cfg.lm_prior_std^2*eye(2*L);
  linit=true(1,L); P.dt=T; P.L=L; P.M=diag([cfg.var_v,cfg.var_w]); P.Q=diag([cfg.var_d,cfg.var_a]);
  x_true=cfg.x0_nom+nz.offset; args=mpc.args; u0=zeros(N,nc); X0=repmat(cfg.x0_nom',N+1,1); k=0;
  while k<maxiter
    x_fb=X(1:3);
    if norm(x_fb(1:2)-cfg.xs(1:2))<cfg.tol, break; end
    args.p=[x_fb;cfg.xs;cfg.safe_buffer];
    args.x0=[reshape(X0',ns*(N+1),1);reshape(u0',nc*N,1)];
    ts=tic;
    sol=mpc.solver('x0',args.x0,'lbx',args.lbx,'ubx',args.ubx,'lbg',args.lbg,'ubg',args.ubg,'p',args.p);
    t_solve(end+1)=toc(ts)*1000;
    solx=full(sol.x); X_sol=reshape(solx(1:ns*(N+1)),ns,N+1)'; u=reshape(solx(ns*(N+1)+1:end),nc,N)'; u_cmd=u(1,:)';
    u_act=u_cmd+nz.u(:,k+1);
    v=u_act(1); w=u_act(2); if abs(w)<1e-9, w=1e-9; end; th=x_true(3);
    x_true=x_true+[v/w*(sin(th+w*T)-sin(th)); v/w*(cos(th)-cos(th+w*T)); w*T]; x_true(3)=wrapToPi(x_true(3));
    z=zeros(2,L); for i=1:L, dx=cfg.lm(1,i)-x_true(1); dy=cfg.lm(2,i)-x_true(2); z(:,i)=[sqrt(dx^2+dy^2);wrapToPi(atan2(dy,dx)-x_true(3))]; end
    z=z+nz.z(:,:,k+1);
    te=tic; [X,Sigma,linit]=mc_ekf_step(X,Sigma,u_cmd,z,linit,P,true); t_ekf(end+1)=toc(te)*1000;
    u0=[u(2:end,:);u(end,:)]; X0=[X_sol(2:end,:);X_sol(end,:)]; k=k+1;
  end
end

T_ms=T*1000; p95=@(v) q_pct(v,95);
fprintf('\n=== Obstacle NMPC per-step timing: %d steps over %d trials, T=%.0f ms ===\n',numel(t_solve),nTrials,T_ms);
fprintf('NMPC solve : mean=%.3f median=%.3f p95=%.3f max=%.3f ms\n',mean(t_solve),median(t_solve),p95(t_solve),max(t_solve));
fprintf('EKF-SLAM   : mean=%.3f median=%.3f p95=%.3f max=%.3f ms\n',mean(t_ekf),median(t_ekf),p95(t_ekf),max(t_ekf));
tot=t_solve(1:numel(t_ekf))+t_ekf;
fprintf('Total/step : mean=%.3f p95=%.3f max=%.3f ms\n',mean(tot),p95(tot),max(tot));
fprintf('Steady-state duty cycle (p95 total / T): %.1f%%\n',100*p95(tot)/T_ms);
fprintf('Worst-case duty cycle  (max total / T): %.1f%%\n',100*max(tot)/T_ms);

function y=q_pct(v,p)
v=sort(v(:)); idx=max(1,ceil(p/100*numel(v))); y=v(idx);
end
