function mpc = mc_build_mpc_dyn(opt)
% MC_BUILD_MPC_DYN  Multiple-shooting NMPC for a single MOVING obstacle.
% The obstacle keep-out is time-varying along the horizon: at node k the robot must
% clear the PREDICTED obstacle centre by a per-node radius (sum of radii + margin).
% Both the predicted centres and the per-node radii are runtime parameters, so the
% same solver serves the static, CV-fixed and covariance-aware variants. A soft
% slack on the obstacle constraint keeps the problem feasible in tight encounters.
%
%   opt.N, opt.T, opt.rob_diam, opt.v_max, opt.omega_max, opt.xy_min/xy_max
%   opt.rho  slack penalty weight (default 1e3)

import casadi.*
if nargin<1 || isempty(opt), opt = struct; end
def = struct('N',14,'T',0.2,'rob_diam',0.3,'v_max',0.6,'omega_max',pi/4, ...
             'xy_min',-2,'xy_max',2,'rho',1e3);
fn = fieldnames(def);
for i=1:numel(fn), if ~isfield(opt,fn{i}), opt.(fn{i})=def.(fn{i}); end, end

N = opt.N; T = opt.T; ns = 3; nc = 2;

x = SX.sym('x'); y = SX.sym('y'); th = SX.sym('theta'); states=[x;y;th];
v = SX.sym('v'); w = SX.sym('omega'); controls=[v;w];
rhs = [v*cos(th); v*sin(th); w];
f = Function('f',{states,controls},{rhs});

U = SX.sym('U',nc,N);
X = SX.sym('X',ns,N+1);
S = SX.sym('S',1,N+1);                       % slack on obstacle constraint per node
% parameters: [x0(3); xref(3); centres(2*(N+1)); radii(N+1)]
P = SX.sym('P', ns + ns + 2*(N+1) + (N+1));
xref = P(4:6);
CB = reshape(P(7 : 6+2*(N+1)), 2, N+1);      % predicted obstacle centres
RB = P(7+2*(N+1) : end);                      % per-node keep-out radii

Wx = diag([1,5,0.1]); Wu = diag([0.5,0.05]);

obj = 0; g = [];
g = [g; X(:,1) - P(1:3)];                     % initial condition
for k = 1:N
    st = X(:,k); con = U(:,k);
    obj = obj + (st-xref)'*Wx*(st-xref) + con'*Wu*con;
    g = [g; X(:,k+1) - (st + T*f(st,con))];   % dynamics (multiple shooting)
end
obj = obj + opt.rho*sumsqr(S);                % slack penalty

% obstacle keep-out (soft): radius_k - ||p_k - centre_k|| - s_k <= 0
for k = 1:N+1
    d = sqrt((X(1,k)-CB(1,k))^2 + (X(2,k)-CB(2,k))^2 + 1e-4);  % eps keeps gradient finite
    g = [g; RB(k) - d - S(k)];
end

OPT = [reshape(X,ns*(N+1),1); reshape(U,nc*N,1); reshape(S,(N+1),1)];
nlp = struct('f',obj,'x',OPT,'g',g,'p',P);
o = struct; o.ipopt.max_iter=150; o.ipopt.print_level=0; o.print_time=0;
o.ipopt.acceptable_tol=1e-8; o.ipopt.acceptable_obj_change_tol=1e-6;
mpc.solver = nlpsol('solver','ipopt',nlp,o);

nX = ns*(N+1); nU = nc*N; nS = N+1;
args = struct;
nEq = ns*(N+1);
args.lbg(1:nEq,1)=0; args.ubg(1:nEq,1)=0;                 % IC + dynamics
args.lbg(nEq+1:nEq+(N+1),1)=-inf; args.ubg(nEq+1:nEq+(N+1),1)=0;   % obstacle

% state bounds
args.lbx(1:3:nX,1)=opt.xy_min; args.ubx(1:3:nX,1)=opt.xy_max;
args.lbx(2:3:nX,1)=opt.xy_min; args.ubx(2:3:nX,1)=opt.xy_max;
args.lbx(3:3:nX,1)=-inf;       args.ubx(3:3:nX,1)=inf;
% control bounds
b = nX;
args.lbx(b+1:2:b+nU,1)=-opt.v_max;     args.ubx(b+1:2:b+nU,1)=opt.v_max;
args.lbx(b+2:2:b+nU,1)=-opt.omega_max; args.ubx(b+2:2:b+nU,1)=opt.omega_max;
% slack bounds (>=0)
c = nX+nU;
args.lbx(c+1:c+nS,1)=0; args.ubx(c+1:c+nS,1)=inf;

mpc.args=args; mpc.f=f; mpc.T=T; mpc.N=N;
mpc.n_states=ns; mpc.n_controls=nc; mpc.n_slack=nS;
mpc.rob_r=opt.rob_diam/2;
end
