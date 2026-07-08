function mpc = mc_build_mpc()
% MC_BUILD_MPC  Build the CasADi single-shooting NMPC once and return a handle
% struct. Same controller as MPC_SLAM_point_stab.m (constrained NMPC, IPOPT),
% factored out so the Monte-Carlo driver can reuse one solver across all trials.

import casadi.*

mpc.T = 0.05;      % sampling time [s]
mpc.N = 10;        % prediction horizon

v_max = 0.6;  v_min = -v_max;
omega_max = pi/4; omega_min = -omega_max;

x = SX.sym('x'); y = SX.sym('y'); theta = SX.sym('theta');
states = [x;y;theta]; n_states = length(states);

v = SX.sym('v'); omega = SX.sym('omega');
controls = [v;omega]; n_controls = length(controls);
rhs = [v*cos(theta); v*sin(theta); omega];

f  = Function('f',{states,controls},{rhs});
U  = SX.sym('U',n_controls,mpc.N);
P  = SX.sym('P',n_states + n_states);    % [current feedback state ; reference]
X  = SX.sym('X',n_states,(mpc.N+1));

% single shooting: roll the state out from the feedback state P(1:3)
X(:,1) = P(1:3);
for k = 1:mpc.N
    st = X(:,k); con = U(:,k);
    X(:,k+1) = st + mpc.T*f(st,con);
end
ff = Function('ff',{U,P},{X});

% cost
Q = diag([1, 5, 0.1]);      % state weights
R = diag([0.5, 0.05]);      % control weights
obj = 0;
for k = 1:mpc.N
    st = X(:,k); con = U(:,k);
    obj = obj + (st-P(4:6))'*Q*(st-P(4:6)) + con'*R*con;
end

% state box constraints (map margins) on x,y over the horizon
g = [];
for k = 1:mpc.N+1
    g = [g; X(1,k)];
    g = [g; X(2,k)];
end

OPT = reshape(U,n_controls*mpc.N,1);
nlp = struct('f',obj,'x',OPT,'g',g,'p',P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level = 0;
opts.print_time = 0;
opts.ipopt.acceptable_tol = 1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

mpc.solver = nlpsol('solver','ipopt',nlp,opts);

args = struct;
args.lbg = -2;  args.ubg = 2;                       % x,y box bounds
args.lbx(1:2:2*mpc.N-1,1) = v_min;     args.lbx(2:2:2*mpc.N,1) = omega_min;
args.ubx(1:2:2*mpc.N-1,1) = v_max;     args.ubx(2:2:2*mpc.N,1) = omega_max;

mpc.args = args;
mpc.ff = ff;
mpc.f  = f;
mpc.n_states = n_states;
mpc.n_controls = n_controls;
mpc.v_max = v_max; mpc.omega_max = omega_max;
end
