function res = mc_run_trial_um2(mode, mpc, cfg, noise)
% MC_RUN_TRIAL_UM2  Waypoint-tour variant of mc_run_trial_um: cfg.wps is a
% 2 x W list of waypoints visited in order (advance when feedback within
% cfg.tol); the LAST waypoint is the goal for the terminal metric. Drift
% suppression is only material over distance, so the tour (~8 m) is what
% separates odometry from SLAM in the unknown-map framing.
%
% UNKNOWN-MAP SLAM framing (start-frame convention):
%   * start pose defines the frame (zero initial pose uncertainty);
%   * landmarks fully unknown, first-observation initialisation,
%     known data association;
%   * waypoints expressed in the start frame.

lm = cfg.lm;  L = size(lm,2);
T  = mpc.T;   N = mpc.N;
wps = cfg.wps;  W = size(wps,2);
wi = 1;                                   % current waypoint index
xs = [wps(:,wi); 0];                      % [x; y; theta_ref]
wp_start_k = 0;                           % step at which current waypoint became active
n  = 3 + 2*L;
maxiter = round(cfg.sim_tim / T);

% --- start-frame initialisation: true = estimate, zero pose uncertainty ---
x_true = cfg.x0_nom;

X = zeros(n,1);  X(1:3) = cfg.x0_nom;
Sigma = zeros(n);
Sigma(4:end,4:end) = cfg.lm_uninf_std^2 * eye(2*L);
linit = false(1,L);

P.dt = T; P.L = L;
P.M = diag([cfg.var_v, cfg.var_w]);
P.Q = diag([cfg.var_d, cfg.var_a]);

args = mpc.args;
u0 = zeros(N, mpc.n_controls);

err_true_ref = zeros(1,maxiter);
err_true_est = zeros(1,maxiter);
sig_t        = zeros(1,maxiter);
solve_t      = zeros(1,maxiter);
traj_true = zeros(3,maxiter+1); traj_true(:,1) = x_true;
traj_est  = zeros(3,maxiter+1); traj_est(:,1)  = X(1:3);

k = 0;
while k < maxiter
    switch mode
        case 'oracle', x_fb = x_true;
        otherwise,     x_fb = X(1:3);
    end
    % arrival test: intermediate waypoints are pass-through (looser tol);
    % only the final goal uses the strict cfg.tol of the terminal metric.
    % A per-waypoint TIMEOUT advances the mission if a waypoint is pursued
    % too long (mission-executive behaviour): an overconfident estimate can
    % hover just outside the tolerance and must not wedge the tour -- the
    % estimation error still pays in the terminal metric.
    if wi == W
        tol_k = cfg.tol;
    else
        tol_k = cfg.tol_wp;
    end
    timed_out = isfield(cfg,'wp_timeout') && wi < W && ...
                (k - wp_start_k)*T > cfg.wp_timeout;
    if norm(x_fb(1:2) - xs(1:2)) < tol_k || timed_out
        if wi == W
            break                          % final waypoint reached
        end
        wi = wi + 1;                       % advance to next leg
        xs = [wps(:,wi); 0];
        wp_start_k = k;
    end

    args.p  = [x_fb; xs];
    args.x0 = reshape(u0', mpc.n_controls*N, 1);
    tsol = tic;
    sol = mpc.solver('x0',args.x0,'lbx',args.lbx,'ubx',args.ubx, ...
                     'lbg',args.lbg,'ubg',args.ubg,'p',args.p);
    solve_t(k+1) = toc(tsol);
    u = reshape(full(sol.x)', mpc.n_controls, N)';
    u_cmd = u(1,:)';

    u_act = u_cmd + noise.u(:,k+1);
    x_true = propagate_exact(x_true, u_act, T);

    if ~strcmp(mode,'oracle')
        z = meas_true(x_true, lm) + noise.z(:,:,k+1);
        % finite sensing range: landmarks are UNKNOWN until first detected
        if isfield(cfg,'r_sense') && cfg.r_sense > 0
            for li = 1:L
                if norm(lm(:,li) - x_true(1:2)) > cfg.r_sense
                    z(:,li) = NaN;
                end
            end
        end
        do_update = strcmp(mode,'slam');
        [X, Sigma, linit] = mc_ekf_step(X, Sigma, u_cmd, z, linit, P, do_update);
    end

    u0 = [u(2:end,:); u(end,:)];

    k = k + 1;
    err_true_ref(k) = norm(x_true(1:2) - wps(:,W));   % vs FINAL goal
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
res.wp_reached    = wi;                                % how far the tour got
res.term_true_ref = norm(x_true(1:2) - wps(:,W));      % vs final waypoint
res.rmse_true_ref = sqrt(mean(err_true_ref.^2));
res.mean_true_est = mean(err_true_est);
res.final_true_est= err_true_est(max(k,1));
res.traj_true     = traj_true(:,1:k+1);
res.traj_est      = traj_est(:,1:k+1);
res.sig_t         = sig_t(1:k);
res.solve_t       = solve_t(1:k);
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
