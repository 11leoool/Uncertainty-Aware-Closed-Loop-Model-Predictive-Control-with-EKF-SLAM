function res = mc_run_trial(mode, mpc, cfg, noise)
% MC_RUN_TRIAL  One closed-loop point-stabilisation trial under a given
% feedback mode, using a PRE-DRAWN noise realisation so that 'oracle',
% 'odom' and 'slam' are compared on the identical scenario.
%
%   mode  : 'oracle' (MPC feedback = true pose)
%           'odom'   (MPC feedback = dead-reckoning estimate, no correction)
%           'slam'   (MPC feedback = EKF-SLAM estimate)
%   mpc   : struct from mc_build_mpc()
%   cfg   : scenario configuration (see run_montecarlo.m)
%   noise : pre-drawn realisation .offset(3x1) .u(2 x maxiter) .z(2 x L x maxiter)
%
% Returns per-trial metrics. The headline metric is term_true_ref: how far the
% REAL robot ends up from the reference (the controller can only act on its
% estimate, so a biased estimate -> biased true endpoint).

lm = cfg.lm;  L = size(lm,2);
T  = mpc.T;   N = mpc.N;
xs = cfg.xs;
n  = 3 + 2*L;
maxiter = round(cfg.sim_tim / T);

% --- true initial pose (offset from nominal) ; estimator starts at nominal ---
x_true = cfg.x0_nom + noise.offset;

X = zeros(n,1);  X(1:3) = cfg.x0_nom;
% Surveyed landmarks: known a priori with an informative prior, refined online.
% This anchors the world frame (so a world-frame goal is well-posed) while the
% filter still estimates & updates the landmark states -> genuine EKF-SLAM.
X(4:end) = cfg.lm(:) + noise.lm(:);
Sigma = zeros(n);
Sigma(1:3,1:3) = diag([cfg.sigma_init_pos^2, cfg.sigma_init_pos^2, cfg.sigma_init_th^2]);
Sigma(4:end,4:end) = cfg.lm_prior_std^2*eye(2*L);
linit = true(1,L);

P.dt = T; P.L = L;
P.M = diag([cfg.var_v, cfg.var_w]);
P.Q = diag([cfg.var_d, cfg.var_a]);

args = mpc.args;
u0 = zeros(N, mpc.n_controls);

err_true_ref = zeros(1,maxiter);
err_true_est = zeros(1,maxiter);
sig_trace    = zeros(1,maxiter);
traj_true = zeros(3,maxiter+1); traj_true(:,1) = x_true;
traj_est  = zeros(3,maxiter+1); traj_est(:,1)  = X(1:3);

k = 0;
while k < maxiter
    switch mode
        case 'oracle', x_fb = x_true;
        otherwise,     x_fb = X(1:3);
    end
    if norm(x_fb(1:2) - xs(1:2)) < cfg.tol
        break
    end

    args.p  = [x_fb; xs];
    args.x0 = reshape(u0', mpc.n_controls*N, 1);
    sol = mpc.solver('x0',args.x0,'lbx',args.lbx,'ubx',args.ubx, ...
                     'lbg',args.lbg,'ubg',args.ubg,'p',args.p);
    u = reshape(full(sol.x)', mpc.n_controls, N)';
    u_cmd = u(1,:)';

    % actuation noise -> applied control drives the TRUE plant
    u_act = u_cmd + noise.u(:,k+1);
    x_true = propagate_exact(x_true, u_act, T);

    % sensor + filter (oracle needs neither)
    if ~strcmp(mode,'oracle')
        z = meas_true(x_true, lm) + noise.z(:,:,k+1);
        do_update = strcmp(mode,'slam');
        [X, Sigma, linit] = mc_ekf_step(X, Sigma, u_cmd, z, linit, P, do_update);
    end

    u0 = [u(2:end,:); u(end,:)];          % warm start

    k = k + 1;
    err_true_ref(k) = norm(x_true(1:2) - xs(1:2));
    if strcmp(mode,'oracle')
        err_true_est(k) = 0;
    else
        err_true_est(k) = norm(x_true(1:2) - X(1:2));
    end
    sig_trace(k) = trace(Sigma(1:3,1:3));
    traj_true(:,k+1) = x_true;
    traj_est(:,k+1)  = X(1:3);
end

err_true_ref = err_true_ref(1:k);
err_true_est = err_true_est(1:k);

res.mode          = mode;
res.steps         = k;
res.term_true_ref = norm(x_true(1:2) - xs(1:2));      % headline metric
res.rmse_true_ref = sqrt(mean(err_true_ref.^2));
res.mean_true_est = mean(err_true_est);
res.final_true_est= err_true_est(max(k,1));
res.traj_true     = traj_true(:,1:k+1);
res.traj_est      = traj_est(:,1:k+1);
res.sig_trace     = sig_trace(1:k);
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
