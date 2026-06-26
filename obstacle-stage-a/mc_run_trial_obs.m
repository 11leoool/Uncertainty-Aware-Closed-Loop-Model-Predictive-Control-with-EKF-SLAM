function res = mc_run_trial_obs(mode, mpc, cfg, noise)
% MC_RUN_TRIAL_OBS  One closed-loop obstacle-avoidance trial (multiple-shooting
% NMPC + EKF-SLAM), same surveyed-landmark / shared-noise setup as
% mc_run_trial.m, plus collision tracking against mpc.obs.
%
%   mode : 'oracle' | 'odom' | 'slam'
%   The MPC plans avoidance using the chosen FEEDBACK pose; collisions are
%   judged on the TRUE pose -> a biased estimate can drive the real robot
%   into the obstacle.
%
% Stage-B hook: if cfg.cov_aware is true and mode=='slam', the obstacle keep-out
% radius is inflated by cfg.gamma*sqrt(lambda_max(Sigma_xy)) (chance constraint).

lm = cfg.lm;  L = size(lm,2);
T  = mpc.T;   N = mpc.N;
xs = cfg.xs;  n = 3 + 2*L;
obs = mpc.obs; K = size(obs,1);
rob_r = mpc.rob_r;
maxiter = round(cfg.sim_tim / T);
cov_aware = isfield(cfg,'cov_aware') && cfg.cov_aware;
gamma = 0; if isfield(cfg,'gamma'), gamma = cfg.gamma; end
base_margin = 0; if isfield(cfg,'safe_buffer'), base_margin = cfg.safe_buffer; end

% --- true initial pose ; estimator starts at nominal ; surveyed landmarks ---
x_true = cfg.x0_nom + noise.offset;
X = zeros(n,1); X(1:3) = cfg.x0_nom; X(4:end) = cfg.lm(:) + noise.lm(:);
Sigma = zeros(n);
Sigma(1:3,1:3) = diag([cfg.sigma_init_pos^2, cfg.sigma_init_pos^2, cfg.sigma_init_th^2]);
Sigma(4:end,4:end) = cfg.lm_prior_std^2*eye(2*L);
linit = true(1,L);

P.dt = T; P.L = L;
P.M = diag([cfg.var_v, cfg.var_w]);
P.Q = diag([cfg.var_d, cfg.var_a]);

args = mpc.args;
u0 = zeros(N, mpc.n_controls);
X0 = repmat((cfg.x0_nom)', N+1, 1);         % multiple-shooting state warm start

err_true_ref = zeros(1,maxiter);
clearance    = inf(1,maxiter);              % true clearance to nearest obstacle
traj_true = zeros(3,maxiter+1); traj_true(:,1) = x_true;
collided = false;

k = 0;
while k < maxiter
    switch mode
        case 'oracle', x_fb = x_true;
        otherwise,     x_fb = X(1:3);
    end
    if norm(x_fb(1:2) - xs(1:2)) < cfg.tol
        break
    end

    % keep-out margin: fixed buffer (all modes) + optional covariance inflation
    % (Stage B, slam only) so the chance constraint grows with pose uncertainty.
    if cov_aware && strcmp(mode,'slam')
        delta = base_margin + gamma * sqrt(max(eig(Sigma(1:2,1:2))));
    else
        delta = base_margin;
    end

    args.p  = [x_fb; xs; delta];
    args.x0 = [reshape(X0', mpc.n_states*(N+1), 1);
               reshape(u0', mpc.n_controls*N, 1)];
    sol = mpc.solver('x0',args.x0,'lbx',args.lbx,'ubx',args.ubx, ...
                     'lbg',args.lbg,'ubg',args.ubg,'p',args.p);

    solx  = full(sol.x);
    X_sol = reshape(solx(1:mpc.n_states*(N+1)), mpc.n_states, N+1)';
    u     = reshape(solx(mpc.n_states*(N+1)+1:end), mpc.n_controls, N)';
    u_cmd = u(1,:)';

    % apply noisy control to the TRUE plant
    u_act = u_cmd + noise.u(:,k+1);
    x_true = propagate_exact(x_true, u_act, T);

    % sensor + filter (oracle uses neither)
    if ~strcmp(mode,'oracle')
        z = meas_true(x_true, lm) + noise.z(:,:,k+1);
        do_update = strcmp(mode,'slam');
        [X, Sigma, linit] = mc_ekf_step(X, Sigma, u_cmd, z, linit, P, do_update);
    end

    % warm starts
    u0 = [u(2:end,:); u(end,:)];
    X0 = [X_sol(2:end,:); X_sol(end,:)];

    k = k + 1;
    err_true_ref(k) = norm(x_true(1:2) - xs(1:2));
    % true clearance to nearest obstacle (negative => physical collision)
    if K > 0
        cl = inf;
        for o = 1:K
            d = norm(x_true(1:2) - obs(o,1:2)') - (rob_r + obs(o,3));
            cl = min(cl, d);
        end
        clearance(k) = cl;
        if cl < 0, collided = true; end
    end
    traj_true(:,k+1) = x_true;
end

err_true_ref = err_true_ref(1:k);
clearance    = clearance(1:k);

res.mode          = mode;
res.steps         = k;
res.term_true_ref = norm(x_true(1:2) - xs(1:2));
res.rmse_true_ref = sqrt(mean(err_true_ref.^2));
res.collided      = collided;
res.min_clear     = min(clearance);
res.traj_true     = traj_true(:,1:k+1);
end

% ---------------------------------------------------------------------------
function x = propagate_exact(x, u, dt)
v = u(1); w = u(2); th = x(3);
if abs(w) < 1e-9, w = 1e-9; end
x = x + [ v/w*( sin(th+w*dt) - sin(th) );
          v/w*( cos(th) - cos(th+w*dt) );
          w*dt ];
x(3) = wrapToPi(x(3));
end

function z = meas_true(x, lm)
L = size(lm,2); z = zeros(2,L);
for i = 1:L
    dx = lm(1,i) - x(1); dy = lm(2,i) - x(2);
    z(:,i) = [ sqrt(dx^2 + dy^2); wrapToPi(atan2(dy,dx) - x(3)) ];
end
end
