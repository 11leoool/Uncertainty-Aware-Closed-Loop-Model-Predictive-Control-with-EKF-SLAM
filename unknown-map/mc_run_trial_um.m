function res = mc_run_trial_um(mode, mpc, cfg, noise)
% MC_RUN_TRIAL_UM  One closed-loop point-stabilisation trial under the
% UNKNOWN-MAP SLAM framing (start-frame convention):
%
%   * The robot's start pose DEFINES the world frame: estimator and true robot
%     both start at cfg.x0_nom with ZERO initial pose uncertainty. (There is no
%     external frame to be wrong about; the global-offset gauge freedom of
%     unknown-map SLAM is thereby fixed.)
%   * Landmarks are FULLY UNKNOWN: initialised from the first observation
%     (inverse measurement model in mc_ekf_step) under an uninformative prior;
%     refined by revisits. Landmarks are assumed distinguishable (known data
%     association).
%   * The goal cfg.xs is expressed in the start frame ("reach the point
%     1.5 m x 1.5 m from where you began").
%   * What SLAM buys over odometry here is DRIFT SUPPRESSION: with actuation
%     noise the dead-reckoning error grows without bound, while re-observing
%     landmarks bounds it.
%
%   mode  : 'oracle' | 'odom' | 'slam'   (same shared-noise discipline)

lm = cfg.lm;  L = size(lm,2);
T  = mpc.T;   N = mpc.N;
xs = cfg.xs;
n  = 3 + 2*L;
maxiter = round(cfg.sim_tim / T);

% --- start-frame initialisation: true = estimate, zero pose uncertainty ---
x_true = cfg.x0_nom;

X = zeros(n,1);  X(1:3) = cfg.x0_nom;
Sigma = zeros(n);
Sigma(4:end,4:end) = cfg.lm_uninf_std^2 * eye(2*L);   % uninformative landmark prior
linit = false(1,L);                                    % landmarks unknown

P.dt = T; P.L = L;
P.M = diag([cfg.var_v, cfg.var_w]);
P.Q = diag([cfg.var_d, cfg.var_a]);

args = mpc.args;
u0 = zeros(N, mpc.n_controls);

err_true_ref = zeros(1,maxiter);
err_true_est = zeros(1,maxiter);
sig_t        = zeros(1,maxiter);      % sqrt(lambda_max(Sigma_xy))
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

    u_act = u_cmd + noise.u(:,k+1);
    x_true = propagate_exact(x_true, u_act, T);

    if ~strcmp(mode,'oracle')
        z = meas_true(x_true, lm) + noise.z(:,:,k+1);
        do_update = strcmp(mode,'slam');
        [X, Sigma, linit] = mc_ekf_step(X, Sigma, u_cmd, z, linit, P, do_update);
    end

    u0 = [u(2:end,:); u(end,:)];

    k = k + 1;
    err_true_ref(k) = norm(x_true(1:2) - xs(1:2));
    if strcmp(mode,'oracle')
        err_true_est(k) = 0;
    else
        err_true_est(k) = norm(x_true(1:2) - X(1:2));
    end
    sig_t(k) = sqrt(max(eig(Sigma(1:2,1:2))));
    traj_true(:,k+1) = x_true;
    traj_est(:,k+1)  = X(1:3);
end

err_true_ref = err_true_ref(1:k);
err_true_est = err_true_est(1:k);

res.mode          = mode;
res.steps         = k;
res.term_true_ref = norm(x_true(1:2) - xs(1:2));
res.rmse_true_ref = sqrt(mean(err_true_ref.^2));
res.mean_true_est = mean(err_true_est);
res.final_true_est= err_true_est(max(k,1));
res.traj_true     = traj_true(:,1:k+1);
res.traj_est      = traj_est(:,1:k+1);
res.sig_t         = sig_t(1:k);
% final landmark-map error (initialised landmarks only; NaN for odom/oracle)
if strcmp(mode,'slam') && any(linit)
    le = [];
    for i = 1:L
        if linit(i)
            le(end+1) = norm(X(2+2*i:3+2*i) - lm(:,i)); %#ok<AGROW>
        end
    end
    res.map_err = mean(le);
else
    res.map_err = NaN;
end
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
