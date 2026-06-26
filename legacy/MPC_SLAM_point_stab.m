% point stabilization + Single shooting
clearvars -except ll
close all
clc

% CasADi v3.4.5
% addpath('C:\Users\mehre\OneDrive\Desktop\CasADi\casadi-windows-matlabR2016a-v3.4.5')
% CasADi v3.5.5
addpath('D:\CODING\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

T = 0.05; % sampling time [s]
N = 10; % prediction horizon
rob_diam = 0.3;

v_max = 0.6; v_min = -v_max;
omega_max = pi/4; omega_min = -omega_max;

x = SX.sym('x'); y = SX.sym('y'); theta = SX.sym('theta');
states = [x;y;theta]; n_states = length(states);

v = SX.sym('v'); omega = SX.sym('omega');
controls = [v;omega]; n_controls = length(controls);
rhs = [v*cos(theta);v*sin(theta);omega]; % system r.h.s

f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U',n_controls,N); % Decision variables (controls)
P = SX.sym('P',n_states + n_states);
% parameters (which include the initial and the reference state of the robot)

X = SX.sym('X',n_states,(N+1));
% A Matrix that represents the states over the optimization problem.

% compute solution symbolically
X(:,1) = P(1:3); % initial state
for k = 1:N
    st = X(:,k);  con = U(:,k);
    f_value  = f(st,con);
    st_next  = st+ (T*f_value);
    X(:,k+1) = st_next;
end
% this function to get the optimal trajectory knowing the optimal solution
ff=Function('ff',{U,P},{X});

obj = 0; % Objective function
g = [];  % constraints vector

Q = zeros(3,3); Q(1,1) = 1;Q(2,2) = 5;Q(3,3) = 0.1; % weighing matrices (states)
R = zeros(2,2); R(1,1) = 0.5; R(2,2) = 0.05; % weighing matrices (controls)
% compute objective
for k=1:N
    st = X(:,k);  con = U(:,k);
    obj = obj+(st-P(4:6))'*Q*(st-P(4:6)) + con'*R*con; % calculate obj
end

% compute constraints
for k = 1:N+1   % box constraints due to the map margins
    g = [g ; X(1,k)];   %state x
    g = [g ; X(2,k)];   %state y
end

% make the decision variables one column vector
OPT_variables = reshape(U,2*N,1);
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);


args = struct;
% inequality constraints (state constraints)
args.lbg = -2;  % lower bound of the states x and y
args.ubg = 2;   % upper bound of the states x and y

% input constraints
args.lbx(1:2:2*N-1,1) = v_min; args.lbx(2:2:2*N,1)   = omega_min;
args.ubx(1:2:2*N-1,1) = v_max; args.ubx(2:2:2*N,1)   = omega_max;


%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SETTING UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
x0 = [0;0; 0.0]; % initial condition.



%true initial condition
x0_true = x0;



xs = [1.5 ; 1.5 ; 0]; % Reference posture.


% Monte Carlo Simulation
%50 times, different starting position N ~ (0,sigma^2)




xx(:,1) = x0; % xx contains the history of states


%initial condition save
x_0_save = x0;

t(1) = t0;

%true state
xx_true(:,1) = x0;
u0 = zeros(N,2);  % two control inputs

sim_tim = 20; % Maximum simulation time

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];
u_cl_SLAM = [];
i = 0;



% load("matlab.mat")
% load("matlab_measurements_noisy_v_omega_d_a.mat")
% the main simulaton loop... it works as long as the error is greater
% than 10^-2 and the number of mpc steps is less than its maximum
% value.

x_location = [-0.5 1 -1.5];
y_location = [0.5 1 0];   % must match the ground-truth landmarks in measure.m / SLAM.m
landmarks = length(x_location);
%initial condition
Sigma = [[zeros(3,3)+0.000001,zeros(3,2*landmarks)];[zeros(2*landmarks,3),1000000000*eye(2*landmarks,2*landmarks)]];
X_p = [0;0;0;0;0;0;0;0;0];
estimate_location = [];
Measurements = [];

error_through_time = [];
error_through_time_truevsestimate = [];
error_through_time_truevsref = [];



% Control-/measurement-noise VARIANCES. These now match the filter model:
%  - control:     M_t = diag([covariance_v, covariance_omega]) in SLAM.m (1e-6)
%  - measurement: Q   = diag([covariance_d, covariance_a])     in SLAM.m (0.01)
covariance_omega = 1e-6;
covariance_v = 1e-6;
covariance_d = 0.01;
covariance_a = 0.01;



[Measurements] = measure([0;0;0],0,Measurements,T);

main_loop = tic;
z(1) = 0;


while(norm((x0-xs),2) > 0.1 && mpciter < sim_tim / T)
    args.p   = [x0;xs]; % set the values of the parameters vector
    args.x0 = reshape(u0',2*N,1); % initial value of the optimization variables
    %tic
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    %toc
    u = reshape(full(sol.x)',2,N)';
    ff_value = ff(u',args.p); % compute OPTIMAL solution TRAJECTORY
    xx1(:,1:3,mpciter+1)= full(ff_value)';

    m = normrnd(0,sqrt(covariance_v));      % normrnd takes std = sqrt(variance)
    n = normrnd(0,sqrt(covariance_omega));

    u_cl= [u_cl ; u(1,:)+[m,n]];
    u_cl_SLAM= [u_cl_SLAM ; u(1,:)];
    t(mpciter+1) = t0;

    t0 = t0+T;
    

    [x_output,y_output,theta_output,Sigma,X_p,estimate_location,Measurements] = SLAM(u_cl,t0,Measurements,Sigma,X_p,estimate_location,T);


    %debug



    x0 = [x_output;y_output;theta_output];


    % true state
    [t0, x0_true, u0_true] = shift_2(T, t0, x0_true,u+[m,n],f);
    
    xx(:,mpciter+2) = x0;

    xx_true(:,mpciter+2) = x0_true;


    error_through_time(mpciter+1) = norm((x0(1:2)-xs(1:2)),2);
    error_through_time_truevsestimate(mpciter+1) = norm((x0_true(1:2)-x0(1:2)),2);
    error_through_time_truevsref(mpciter+1) = norm((x0_true(1:2)-xs(1:2)),2);
    mpciter
    mpciter = mpciter + 1;
    

end;
main_loop_time = toc(main_loop)

ss_error = norm((x0(1:2)-xs(1:2)),2);
ss_error_truevsestimate = norm((x0_true(1:2)-x0(1:2)),2);
ss_error_truevsref = norm((x0_true(1:2)-xs(1:2)),2);
ss = [ss_error;ss_error_truevsestimate;ss_error_truevsref];

ss_error_angle = norm((x0(3)-xs(3)),2);
ss_error_truevsestimate_angle = norm((x0_true(3)-x0(3)),2);
ss_error_truevsref_angle = norm((x0_true(3)-xs(3)),2);
ss_angle = [ss_error_angle;ss_error_truevsestimate_angle;ss_error_truevsref_angle];

% Localisation error of the EKF-SLAM estimate vs the true state.
% (error_through_time_truevsestimate stores per-step Euclidean norms.)
MAE  = mean(error_through_time_truevsestimate);                 % mean absolute error
RMSE = sqrt(mean(error_through_time_truevsestimate.^2));        % root-mean-square error


%KF with true state
%RF with true state
%RF with KF estiamte
%repeat the simulation for 50 times take average error
%Monte Carlo Simulation
%random starting point,measurement/motion noise
%true landmarks and estimated landmarks

MPC_SLAM_stab_draw (t,xx,xx1,u_cl,xs,N,rob_diam,xx_true,u_cl_SLAM) % a drawing function
hold off

figure(999)
plot(1:1:mpciter,error_through_time);
hold on
grid on
plot(1:1:mpciter,error_through_time_truevsestimate);
plot(1:1:mpciter,error_through_time_truevsref);
legend("estimated position vs reference position","estimated position vs true position","true position vs reference position");
title("distance Error vs time");
xlabel("step time(0.2s)");
ylabel("error (meters)");



% saveas(gcf,'SLAM.fig')