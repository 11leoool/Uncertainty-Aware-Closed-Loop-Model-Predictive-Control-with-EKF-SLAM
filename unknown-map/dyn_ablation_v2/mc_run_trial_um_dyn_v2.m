function res = mc_run_trial_um_dyn_v2(mode, mpc, cfg, noise)
% MC_RUN_TRIAL_UM_DYN  Stage 4 of the unknown-map refactor: patrol tour with
% a MOVING obstacle, everything unknown until detected (finite sensing range):
% landmarks enter the map on first detection; the obstacle track is created on
% first detection and coasts (prediction-only, growing covariance) when out of
% range. Robot pose from EKF-SLAM (start-frame convention).
%
%   mode :
%     'oracle'   true pose + clairvoyant true obstacle future, fixed margin
%     'static'   SLAM pose + obstacle frozen at current estimate, fixed margin
%     'cv_fixed' SLAM pose + CV-predicted track, fixed margin
%     'cv_cov'   SLAM pose + CV-predicted track, covariance-aware margin
%                delta_k = d0 + gamma*sqrt(lambda_max(P_xy + Sigma_obs,k))
%
% Collisions judged against the TRUE moving obstacle. Constraint inactive
% until the obstacle is first detected (except oracle, which is clairvoyant).

lm = cfg.lm;  L = size(lm,2);
T  = mpc.T;   N = mpc.N;
wps = cfg.wps;  W = size(wps,2);
wi = 1;  xs = [wps(:,wi); 0];  wp_start_k = 0;
n_aug = 3 + 2*L;
rob_r = mpc.rob_r;  obs_r = cfg.obs_r;
d0 = cfg.safe_buffer;  gamma = cfg.gamma;
% cfg.fixed_extra: CONSTANT margin inflation for the matched-mean ablation
% (applies to the cv_fixed strategy in place of the covariance term)
fixed_extra = 0; if isfield(cfg,'fixed_extra'), fixed_extra = cfg.fixed_extra; end
maxiter = round(cfg.sim_tim / T);

% --- TRUE obstacle trajectory (CV + accel noise), with horizon look-ahead ---
F4 = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];
G4 = [T^2/2 0; 0 T^2/2; T 0; 0 T];
ntot = maxiter + N + 2;
Otrue = zeros(2, ntot);
ost = [cfg.o0; cfg.vo0];
for nn = 1:ntot
    Otrue(:,nn) = ost(1:2);
    ost = F4*ost + G4*noise.oacc(:,nn);
end

% --- start-frame robot init, unknown landmarks ---
x_true = cfg.x0_nom;
X = zeros(n_aug,1); X(1:3) = cfg.x0_nom;
Sig = zeros(n_aug);
Sig(4:end,4:end) = cfg.lm_uninf_std^2*eye(2*L);
linit = false(1,L);
Pp.dt=T; Pp.L=L; Pp.M=diag([cfg.var_v,cfg.var_w]); Pp.Q=diag([cfg.var_d,cfg.var_a]);

% --- obstacle tracker (uninitialised until first detection) ---
% Track management: after cfg.track_drop seconds without a detection the
% track is DROPPED (constraint deactivates until re-detection). Without
% this, the coasting covariance grows unboundedly and the covariance-aware
% margin inflates a stale phantom bubble that can paralyse the mission.
o_est = zeros(4,1); So = eye(4); oinit = false;
miss_cnt = 0;
drop_steps = inf;
if isfield(cfg,'track_drop'), drop_steps = round(cfg.track_drop / T); end
Pt.dt=T; Pt.sa2=cfg.obs_sa2_filter; Pt.Ro=diag([cfg.var_d,cfg.var_a]);
Qo = G4*(Pt.sa2*eye(2))*G4';

ns=mpc.n_states; nc=mpc.n_controls; nX=ns*(N+1); nU=nc*N; nS=N+1;
args = mpc.args;
X0 = repmat(cfg.x0_nom,1,N+1); u0 = zeros(nc,N); S0 = zeros(1,N+1);

clr = inf(1,maxiter);
infl_t = nan(1,maxiter);              % node-0 inflation; NaN when constraint inactive
active_t = false(1,maxiter);          % was an obstacle constraint active at MPC build
slack_t = nan(1,maxiter);             % max soft-constraint slack per step
infl_now = NaN; act_now = false;      % per-step status, set in the horizon block
sig_t  = zeros(1,maxiter);            % reported robot pose sigma per step
traj = zeros(3,maxiter+1); traj(:,1)=x_true;
traj_est = zeros(3,maxiter+1); traj_est(:,1)=X(1:3);
otraj = zeros(2,maxiter+1); otraj(:,1)=Otrue(:,1);
collided = false; detect_step = NaN; k = 0;

while k < maxiter
    nn = k+1;
    if strcmp(mode,'oracle'), x_fb = x_true; else, x_fb = X(1:3); end

    if wi == W, tol_k = cfg.tol; else, tol_k = cfg.tol_wp; end
    timed_out = isfield(cfg,'wp_timeout') && wi < W && ...
                (k - wp_start_k)*T > cfg.wp_timeout;
    if norm(x_fb(1:2) - xs(1:2)) < tol_k || timed_out
        if wi == W, break; end
        wi = wi + 1;  xs = [wps(:,wi); 0];  wp_start_k = k;
    end

    % ---- predicted obstacle horizon ----
    CB = zeros(2,N+1); RB = zeros(N+1,1);
    if strcmp(mode,'oracle')
        for j = 0:N, CB(:,j+1) = Otrue(:, nn+j); RB(j+1) = rob_r+obs_r+d0; end
        infl_now = 0; act_now = true;
    elseif ~oinit
        CB(:) = 50; RB(:) = 1e-3;          % not yet detected: inactive
        infl_now = NaN; act_now = false;
    elseif strcmp(mode,'static')
        for j = 0:N, CB(:,j+1) = o_est(1:2); RB(j+1) = rob_r+obs_r+d0; end
        infl_now = 0; act_now = true;
    else  % cv_fixed / cv_cov
        Pxy = Sig(1:2,1:2);
        oj = o_est; Sj = So;
        for j = 0:N
            CB(:,j+1) = oj(1:2);
            if strcmp(mode,'cv_cov')
                dlt = d0 + gamma*sqrt(max(eig(Pxy + Sj(1:2,1:2))));
            else
                dlt = d0 + fixed_extra;
            end
            RB(j+1) = rob_r+obs_r+dlt;
            if j == 0, infl_now = dlt - d0; act_now = true; end
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

    % ---- true plant ----
    u_act = u_cmd + noise.u(:,nn);
    x_true = propagate_exact(x_true, u_act, T);

    % ---- sense + filters (oracle uses neither) ----
    if ~strcmp(mode,'oracle')
        z_lm = meas_lm(x_true, lm) + noise.z(:,:,nn);
        for li = 1:L
            if norm(lm(:,li) - x_true(1:2)) > cfg.r_sense
                z_lm(:,li) = NaN;
            end
        end
        [X, Sig, linit] = mc_ekf_step(X, Sig, u_cmd, z_lm, linit, Pp, true);

        % obstacle detection only within sensing range
        if norm(Otrue(:,nn+1) - x_true(1:2)) <= cfg.r_sense
            z_o = meas_pt(x_true, Otrue(:,nn+1)) + noise.oz(:,nn);
            if ~oinit && isnan(detect_step), detect_step = k+1; end
            miss_cnt = 0;
        else
            z_o = [];
            miss_cnt = miss_cnt + 1;
        end
        if oinit && miss_cnt >= drop_steps
            oinit = false;                 % stale track: drop until re-detected
            o_est = zeros(4,1); So = eye(4);
        end
        [o_est, So, oinit] = cv_tracker_um(o_est, So, z_o, X(1:3), Sig(1:3,1:3), oinit, Pt);
    end

    % warm starts
    X0 = [X_sol(:,2:end), X_sol(:,end)];
    u0 = [u(:,2:end), u(:,end)];
    S0 = zeros(1,N+1);

    k = k+1;
    clr(k) = norm(x_true(1:2)-Otrue(:,nn+1)) - (rob_r+obs_r);
    if clr(k) < 0, collided = true; end
    active_t(k) = act_now;
    if act_now, infl_t(k) = infl_now; end    % stays NaN when inactive
    try, slack_t(k) = max(solx(nX+nU+1:nX+nU+nS)); catch, end
    sig_t(k) = sqrt(max(eig(Sig(1:2,1:2))));    % reported pose sigma
    traj(:,k+1) = x_true;  traj_est(:,k+1) = X(1:3);
    otraj(:,k+1) = Otrue(:,nn+1);
end

res.mode = mode; res.steps = k;
res.wp_reached = wi;
res.term_true_ref = norm(x_true(1:2)-wps(:,W));
res.collided = collided;
res.min_clear = min(clr(1:max(k,1)));
res.detect_step = detect_step;
res.path_len = sum(sqrt(sum(diff(traj(1:2,1:k+1),1,2).^2,1)));
res.traj = traj(:,1:k+1); res.traj_est = traj_est(:,1:k+1);
res.otraj = otraj(:,1:k+1);
res.mean_loc_err = mean(vecnorm(traj(1:2,1:k+1) - traj_est(1:2,1:k+1)));
ia = infl_t(1:k); aa = active_t(1:k);
assert(all(isnan(ia(~aa))), 'inflation recorded while constraint inactive');
res.active_frac = mean(aa);
if any(aa)
    res.mean_inflation_act = mean(ia(aa)); res.peak_inflation = max(ia(aa));
    res.med_inflation_act = median(ia(aa)); res.iqr_inflation_act = iqr(ia(aa));
else
    res.mean_inflation_act = 0; res.peak_inflation = 0;
    res.med_inflation_act = 0; res.iqr_inflation_act = 0;
end
res.mean_inflation_all = sum(ia(aa))/max(k,1);
assert(res.mean_inflation_all <= res.active_frac*res.peak_inflation + 1e-12, ...
    'stale-carry invariant violated: mean_all > active_frac*peak');
res.slack_freq = mean(slack_t(1:k) > 1e-6); res.slack_max = max([slack_t(1:k), 0]);
res.mean_inflation = res.mean_inflation_all;   % legacy name, corrected definition
res.mean_sigma = mean(sig_t(1:max(k,1)));
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
