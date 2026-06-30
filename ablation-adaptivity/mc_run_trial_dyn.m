function res = mc_run_trial_dyn(mode, mpc, cfg, noise)
% MC_RUN_TRIAL_DYN  One closed-loop trial with a MOVING obstacle.
% Robot pose is estimated by EKF-SLAM; the obstacle is tracked by a CV-EKF
% (cv_tracker_step). The MPC avoids the PREDICTED obstacle track over the horizon.
%
%   mode :
%     'oracle'   true robot pose + clairvoyant true obstacle future, fixed margin
%     'static'   SLAM pose + obstacle frozen at current estimate, fixed margin
%     'cv_fixed' SLAM pose + CV-predicted obstacle track, fixed margin
%     'cv_cov'   SLAM pose + CV-predicted track, covariance-aware margin
%                delta_k = delta0 + gamma*sqrt(lambda_max(P_xy + Sigma_obs,k))
%
% Collisions are judged against the TRUE moving obstacle.

lm = cfg.lm;  L = size(lm,2);
T  = mpc.T;   N = mpc.N;
xs = cfg.xs;  n_aug = 3 + 2*L;
rob_r = mpc.rob_r;  obs_r = cfg.obs_r;
d0 = cfg.safe_buffer;  gamma = cfg.gamma;
fixed_extra = 0; if isfield(cfg,'fixed_extra'), fixed_extra = cfg.fixed_extra; end
infl_sum = 0; infl_n = 0;          % accumulate applied inflation (delta_k - d0) over nodes/steps
maxiter = round(cfg.sim_tim / T);

% --- pre-compute the TRUE obstacle trajectory (CV + accel noise), with look-ahead
F4 = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];
G4 = [T^2/2 0; 0 T^2/2; T 0; 0 T];
ntot = maxiter + N + 2;
Otrue = zeros(2, ntot);
ost = [cfg.o0; cfg.vo0];
for nn = 1:ntot
    Otrue(:,nn) = ost(1:2);
    ost = F4*ost + G4*noise.oacc(:,nn);     % true accel perturbation
end

% --- robot truth + EKF-SLAM init (surveyed landmarks) ---
x_true = cfg.x0_nom + noise.offset;
X = zeros(n_aug,1); X(1:3) = cfg.x0_nom; X(4:end) = cfg.lm(:) + noise.lm(:);
Sig = zeros(n_aug);
Sig(1:3,1:3) = diag([cfg.sigma_init_pos^2, cfg.sigma_init_pos^2, cfg.sigma_init_th^2]);
Sig(4:end,4:end) = cfg.lm_prior_std^2*eye(2*L);
linit = true(1,L);
Pp.dt=T; Pp.L=L; Pp.M=diag([cfg.var_v,cfg.var_w]); Pp.Q=diag([cfg.var_d,cfg.var_a]);

% --- obstacle tracker init ---
o_est = zeros(4,1); So = eye(4); oinit = false;
Pt.dt=T; Pt.sa2=cfg.obs_sa2_filter; Pt.Ro=diag([cfg.var_d,cfg.var_a]);
Qo = G4*(Pt.sa2*eye(2))*G4';

ns=mpc.n_states; nc=mpc.n_controls; nX=ns*(N+1); nU=nc*N; nS=N+1;
args = mpc.args;
X0 = repmat(cfg.x0_nom,1,N+1); u0 = zeros(nc,N); S0 = zeros(1,N+1);

err_tr = zeros(1,maxiter); clr = inf(1,maxiter);
traj = zeros(3,maxiter+1); traj(:,1)=x_true;
otraj = zeros(2,maxiter+1); otraj(:,1)=Otrue(:,1);
collided = false; k = 0;

while k < maxiter
    nn = k+1;
    if strcmp(mode,'oracle'), x_fb = x_true; else, x_fb = X(1:3); end
    if norm(x_fb(1:2)-xs(1:2)) < cfg.tol, break; end

    % ---- build predicted obstacle horizon (centres CB, per-node radius RB) ----
    CB = zeros(2,N+1); RB = zeros(N+1,1);
    if strcmp(mode,'oracle')
        for j = 0:N, CB(:,j+1) = Otrue(:, nn+j); RB(j+1) = rob_r+obs_r+d0; end
    elseif strcmp(mode,'static')
        for j = 0:N, CB(:,j+1) = o_est(1:2); RB(j+1) = rob_r+obs_r+d0; end
    else  % cv_fixed / cv_cov : CV roll-out of the tracker
        Pxy = Sig(1:2,1:2);
        oj = o_est; Sj = So;
        for j = 0:N
            CB(:,j+1) = oj(1:2);
            if strcmp(mode,'cv_cov')
                Sxy = Sj(1:2,1:2);
                dlt = d0 + gamma*sqrt(max(eig(Pxy + Sxy)));
            elseif strcmp(mode,'cv_fixedmatch')
                dlt = d0 + fixed_extra;            % constant margin matched to cv_cov mean
            else
                dlt = d0;
            end
            RB(j+1) = rob_r+obs_r+dlt;
            infl_sum = infl_sum + dlt - d0; infl_n = infl_n + 1;
            oj = F4*oj; Sj = F4*Sj*F4' + Qo;
        end
    end

    % ---- solve MPC ----
    args.p = [x_fb; xs; CB(:); RB];
    args.x0 = [reshape(X0,nX,1); reshape(u0,nU,1); reshape(S0,nS,1)];
    sol = mpc.solver('x0',args.x0,'lbx',args.lbx,'ubx',args.ubx, ...
                     'lbg',args.lbg,'ubg',args.ubg,'p',args.p);
    solx = full(sol.x);
    X_sol = reshape(solx(1:nX), ns, N+1);
    u     = reshape(solx(nX+1:nX+nU), nc, N);
    u_cmd = u(:,1);

    % ---- apply control to the true plant ----
    u_act = u_cmd + noise.u(:,nn);
    x_true = propagate_exact(x_true, u_act, T);

    % ---- sense + update filters (oracle uses neither) ----
    if ~strcmp(mode,'oracle')
        z_lm = meas_lm(x_true, lm) + noise.z(:,:,nn);
        [X, Sig, linit] = mc_ekf_step(X, Sig, u_cmd, z_lm, linit, Pp, true);
        z_o = meas_pt(x_true, Otrue(:,nn+1)) + noise.oz(:,nn);
        [o_est, So, oinit] = cv_tracker_step(o_est, So, z_o, X(1:3), oinit, Pt);
    end

    % warm starts
    X0 = [X_sol(:,2:end), X_sol(:,end)];
    u0 = [u(:,2:end), u(:,end)];
    S0 = zeros(1,N+1);

    k = k+1;
    err_tr(k) = norm(x_true(1:2)-xs(1:2));
    clr(k) = norm(x_true(1:2)-Otrue(:,nn+1)) - (rob_r+obs_r);
    if clr(k) < 0, collided = true; end
    traj(:,k+1) = x_true;  otraj(:,k+1) = Otrue(:,nn+1);
end

err_tr = err_tr(1:k); clr = clr(1:k);
res.mode = mode; res.steps = k;
res.term_true_ref = norm(x_true(1:2)-xs(1:2));
res.rmse_true_ref = sqrt(mean(err_tr.^2));
res.collided = collided;
res.min_clear = min(clr);
res.path_len = sum(sqrt(sum(diff(traj(1:2,1:k+1),1,2).^2,1)));
res.infl_mean = infl_sum / max(infl_n,1);
res.traj = traj(:,1:k+1); res.otraj = otraj(:,1:k+1);
end

% ---------------------------------------------------------------------------
function x = propagate_exact(x,u,dt)
v=u(1); w=u(2); th=x(3);
if abs(w)<1e-9, w=1e-9; end
x = x + [ v/w*(sin(th+w*dt)-sin(th)); v/w*(cos(th)-cos(th+w*dt)); w*dt ];
x(3)=wrapToPi(x(3));
end

function z = meas_lm(x,lm)
L=size(lm,2); z=zeros(2,L);
for i=1:L
    dx=lm(1,i)-x(1); dy=lm(2,i)-x(2);
    z(:,i)=[sqrt(dx^2+dy^2); wrapToPi(atan2(dy,dx)-x(3))];
end
end

function z = meas_pt(x,p)
dx=p(1)-x(1); dy=p(2)-x(2);
z=[sqrt(dx^2+dy^2); wrapToPi(atan2(dy,dx)-x(3))];
end
