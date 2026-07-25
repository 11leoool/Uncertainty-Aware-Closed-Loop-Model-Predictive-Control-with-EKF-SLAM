function res = mc_run_trial_um_obs(mode, mpc, cfg, noise)
% MC_RUN_TRIAL_UM_OBS  Stage 2 of the unknown-map refactor: waypoint tour
% with a STATIC obstacle that is UNKNOWN UNTIL DETECTED (finite sensing
% range), estimated online by static_obs_ekf from the robot's own pose
% estimate. Landmarks likewise unknown until within cfg.r_sense.
%
%   mode : 'oracle'     true pose feedback, clairvoyant true obstacle
%          'odom_fixed' dead-reckoning pose, sensed obstacle, fixed margin
%          'slam_fixed' EKF-SLAM pose,      sensed obstacle, fixed margin
%          'slam_cov'   EKF-SLAM pose,      sensed obstacle, covariance-aware
%                       margin  delta = d0 + gamma*sqrt(lmax(P_xy + P_obs))
%
% Before the obstacle is first detected the keep-out constraint is inactive
% (dummy far-away centre). Collisions are judged on TRUE robot vs TRUE
% obstacle geometry.

lm = cfg.lm;  L = size(lm,2);
T  = mpc.T;   N = mpc.N;
wps = cfg.wps;  W = size(wps,2);
wi = 1;  xs = [wps(:,wi); 0];
wp_start_k = 0;
n  = 3 + 2*L;
maxiter = round(cfg.sim_tim / T);
ns = mpc.n_states; nc = mpc.n_controls;
nX = ns*(N+1); nU = nc*N; nS = mpc.n_slack;

o_true = cfg.obs(1:2)';  r_obs = cfg.obs(3);
rob_r  = mpc.rob_r;
gamma  = cfg.gamma;  d0 = cfg.safe_buffer;
% cfg.filt_scale: filter-consistency knob. The filter ASSUMES noise variances
% scaled by filt_scale while the TRUE noise draws (in the runner) are unchanged.
% filt_scale < 1 -> overconfident filter (reports P smaller than the true error);
% filt_scale > 1 -> conservative filter; filt_scale = 1 -> consistent (default).
fs = 1; if isfield(cfg,'filt_scale'), fs = cfg.filt_scale; end
Rz = diag([cfg.var_d, cfg.var_a]) * fs;

slam_fb  = any(strcmp(mode, {'slam_fixed','slam_cov'}));
use_est_obs = ~strcmp(mode,'oracle');

% --- start-frame initialisation ---
x_true = cfg.x0_nom;
X = zeros(n,1);  X(1:3) = cfg.x0_nom;
Sigma = zeros(n);
Sigma(4:end,4:end) = cfg.lm_uninf_std^2 * eye(2*L);
linit = false(1,L);

P.dt = T; P.L = L;
P.M = diag([cfg.var_v, cfg.var_w]) * fs;   % filter's ASSUMED process noise
P.Q = diag([cfg.var_d, cfg.var_a]) * fs;   % filter's ASSUMED measurement noise

o_hat = [];  Po = [];  detected = false;   % obstacle estimator state

args = mpc.args;
X0 = repmat(cfg.x0_nom,1,N+1); u0 = zeros(nc,N); S0 = zeros(1,N+1);

clearance = inf(1,maxiter);
solve_t = zeros(1,maxiter);
sig_t   = zeros(1,maxiter);
errloc_t= zeros(1,maxiter);
delta_t = zeros(1,maxiter);
obserr_t= nan(1,maxiter);
traj_true = zeros(3,maxiter+1); traj_true(:,1) = x_true;
collided = false; detect_step = NaN;

k = 0;
while k < maxiter
    switch mode
        case 'oracle', x_fb = x_true;
        otherwise,     x_fb = X(1:3);
    end
    if wi == W, tol_k = cfg.tol; else, tol_k = cfg.tol_wp; end
    timed_out = isfield(cfg,'wp_timeout') && wi < W && ...
                (k - wp_start_k)*T > cfg.wp_timeout;
    if norm(x_fb(1:2) - xs(1:2)) < tol_k || timed_out
        if wi == W, break; end
        wi = wi + 1;  xs = [wps(:,wi); 0];
        wp_start_k = k;
    end

    % ----- keep-out margin -----
    % (cfg.fixed_extra: CONSTANT inflation for the matched-fixed ablation)
    if strcmp(mode,'slam_cov') && detected
        delta = d0 + gamma*sqrt(max(eig(Sigma(1:2,1:2) + Po)));
    elseif isfield(cfg,'fixed_extra') && cfg.fixed_extra > 0 && ~strcmp(mode,'oracle')
        delta = d0 + cfg.fixed_extra;
    else
        delta = d0;
    end

    % ----- obstacle constraint parameters -----
    if strcmp(mode,'oracle')
        cb = repmat(o_true,1,N+1);
        rb = (rob_r + r_obs + delta)*ones(N+1,1);
    elseif detected
        cb = repmat(o_hat,1,N+1);
        rb = (rob_r + r_obs + delta)*ones(N+1,1);
    else
        cb = repmat([50;50],1,N+1);        % not yet detected: inactive
        rb = 1e-3*ones(N+1,1);
    end

    args.p  = [x_fb; xs; reshape(cb,2*(N+1),1); rb];
    args.x0 = [reshape(X0,nX,1); reshape(u0,nU,1); reshape(S0,nS,1)];
    tsol = tic;
    sol = mpc.solver('x0',args.x0,'lbx',args.lbx,'ubx',args.ubx, ...
                     'lbg',args.lbg,'ubg',args.ubg,'p',args.p);
    solve_t(k+1) = toc(tsol);
    solx  = full(sol.x);
    X_sol = reshape(solx(1:nX), ns, N+1);
    u     = reshape(solx(nX+1:nX+nU), nc, N);
    u_cmd = u(:,1);

    % ----- true plant -----
    u_act = u_cmd + noise.u(:,k+1);
    x_true = propagate_exact(x_true, u_act, T);

    % ----- landmark sensing + pose filter -----
    if ~strcmp(mode,'oracle')
        z = meas_true(x_true, lm) + noise.z(:,:,k+1);
        for li = 1:L
            if norm(lm(:,li) - x_true(1:2)) > cfg.r_sense
                z(:,li) = NaN;
            end
        end
        do_update = slam_fb;
        [X, Sigma, linit] = mc_ekf_step(X, Sigma, u_cmd, z, linit, P, do_update);
    end

    % ----- obstacle sensing + estimator (all non-oracle modes) -----
    if use_est_obs && norm(o_true - x_true(1:2)) <= cfg.r_sense
        dxo = o_true(1)-x_true(1); dyo = o_true(2)-x_true(2);
        zo = [sqrt(dxo^2+dyo^2); wrapToPi(atan2(dyo,dxo)-x_true(3))] + noise.zo(:,k+1);
        if ~detected
            [o_hat, Po] = static_obs_ekf('init', [], [], zo, X(1:3), Sigma(1:3,1:3), Rz);
            detected = true;  detect_step = k+1;
        else
            [o_hat, Po] = static_obs_ekf('update', o_hat, Po, zo, X(1:3), Sigma(1:3,1:3), Rz);
        end
    end

    % ----- warm starts -----
    u0 = [u(:,2:end), u(:,end)];
    X0 = [X_sol(:,2:end), X_sol(:,end)];
    S0 = zeros(1,N+1);

    k = k + 1;
    sig_t(k)   = sqrt(max(eig(Sigma(1:2,1:2))));
    errloc_t(k)= norm(X(1:2) - x_true(1:2));    % TRUE localization error
    delta_t(k) = delta - d0;
    if detected, obserr_t(k) = norm(o_hat - o_true); end
    cl = norm(x_true(1:2) - o_true) - (rob_r + r_obs);
    clearance(k) = cl;
    if cl < 0, collided = true; end
    traj_true(:,k+1) = x_true;
end

res.mode          = mode;
res.steps         = k;
res.wp_reached    = wi;
res.term_true_ref = norm(x_true(1:2) - wps(:,W));
res.collided      = collided;
res.min_clear     = min(clearance(1:max(k,1)));
res.detect_step   = detect_step;
res.obs_err_final = obserr_t(max(k,1));
res.traj_true     = traj_true(:,1:k+1);
res.sig_t         = sig_t(1:k);
res.solve_t       = solve_t(1:k);
res.delta_t       = delta_t(1:k);
res.obserr_t      = obserr_t(1:k);
res.mean_inflation= mean(delta_t(1:k));
res.mean_loc_err  = mean(errloc_t(1:k));
res.mean_sigma    = mean(sig_t(1:k));
res.path_len      = sum(vecnorm(diff(traj_true(1:2,1:k+1),1,2)));
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
