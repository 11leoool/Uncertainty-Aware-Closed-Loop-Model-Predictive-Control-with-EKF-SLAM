function mpc = mc_build_mpc_obs(opt)
% MC_BUILD_MPC_OBS  Multiple-shooting NMPC with disk obstacle-avoidance
% constraints. Used by the obstacle Monte-Carlo study. The single safety
% margin P(7) (=delta) is added to every obstacle keep-out radius at runtime,
% so the covariance-aware variant (Stage B) only needs to pass delta>0.
%
%   opt.N, opt.T          horizon / sampling
%   opt.rob_diam          robot diameter
%   opt.v_max, opt.omega_max
%   opt.obs (K x 3)       obstacle disks, each row [ox oy r_obs]
%   opt.xy_min/xy_max     map box bounds

import casadi.*
if nargin < 1 || isempty(opt), opt = struct; end
def = struct('N',14,'T',0.2,'rob_diam',0.3,'v_max',0.6,'omega_max',pi/4, ...
             'obs',zeros(0,3),'xy_min',-2,'xy_max',2);
fn = fieldnames(def);
for i = 1:numel(fn), if ~isfield(opt,fn{i}), opt.(fn{i}) = def.(fn{i}); end, end

N = opt.N; T = opt.T; rob_r = opt.rob_diam/2;
obs = opt.obs; K = size(obs,1);
ns = 3; nc = 2;

x = SX.sym('x'); y = SX.sym('y'); th = SX.sym('theta'); states = [x;y;th];
v = SX.sym('v'); w = SX.sym('omega'); controls = [v;w];
rhs = [v*cos(th); v*sin(th); w];
f = Function('f',{states,controls},{rhs});

U = SX.sym('U',nc,N);
X = SX.sym('X',ns,N+1);
P = SX.sym('P', ns + ns + 1);          % [feedback state ; reference ; safety margin delta]

Q = diag([1,5,0.1]); R = diag([0.5,0.05]);

obj = 0; g = [];
g = [g; X(:,1) - P(1:3)];               % initial-condition constraint
for k = 1:N
    st = X(:,k); con = U(:,k);
    obj = obj + (st-P(4:6))'*Q*(st-P(4:6)) + con'*R*con;
    st_next = X(:,k+1);
    g = [g; st_next - (st + T*f(st,con))];   % multiple-shooting dynamics
end

% obstacle keep-out: (rob_r + obs_r + delta) - dist <= 0  ->  dist >= safe+delta
for o = 1:K
    ox = obs(o,1); oy = obs(o,2); orad = obs(o,3);
    for k = 1:N+1
        d = sqrt((X(1,k)-ox)^2 + (X(2,k)-oy)^2);
        g = [g; (rob_r + orad + P(7)) - d];
    end
end

OPT = [reshape(X,ns*(N+1),1); reshape(U,nc*N,1)];
nlp = struct('f',obj,'x',OPT,'g',g,'p',P);

o = struct;
o.ipopt.max_iter = 100; o.ipopt.print_level = 0; o.print_time = 0;
o.ipopt.acceptable_tol = 1e-8; o.ipopt.acceptable_obj_change_tol = 1e-6;
mpc.solver = nlpsol('solver','ipopt',nlp,o);

args = struct;
nEq = ns*(N+1);
args.lbg(1:nEq,1) = 0;  args.ubg(1:nEq,1) = 0;             % IC + dynamics (equality)
nIneq = (N+1)*K;
if nIneq > 0
    args.lbg(nEq+1:nEq+nIneq,1) = -inf;                    % obstacle (inequality)
    args.ubg(nEq+1:nEq+nIneq,1) = 0;
end

% state bounds (x,y boxed; theta free)
args.lbx(1:3:ns*(N+1),1) = opt.xy_min; args.ubx(1:3:ns*(N+1),1) = opt.xy_max;
args.lbx(2:3:ns*(N+1),1) = opt.xy_min; args.ubx(2:3:ns*(N+1),1) = opt.xy_max;
args.lbx(3:3:ns*(N+1),1) = -inf;       args.ubx(3:3:ns*(N+1),1) = inf;
% control bounds
b = ns*(N+1);
args.lbx(b+1:2:b+nc*N,1) = -opt.v_max;     args.ubx(b+1:2:b+nc*N,1) = opt.v_max;
args.lbx(b+2:2:b+nc*N,1) = -opt.omega_max; args.ubx(b+2:2:b+nc*N,1) = opt.omega_max;

mpc.args = args; mpc.f = f; mpc.T = T; mpc.N = N;
mpc.n_states = ns; mpc.n_controls = nc;
mpc.obs = obs; mpc.rob_diam = opt.rob_diam; mpc.rob_r = rob_r;
end
