function [X, Sigma, linit] = mc_ekf_step(X, Sigma, u_cmd, z, linit, P, do_update)
% MC_EKF_STEP  One EKF-SLAM predict (+ optional sequential measurement update).
%
%   X       (3+2L x 1)  state estimate [x;y;theta; l1x;l1y; ...]
%   Sigma   (n x n)     covariance
%   u_cmd   (2 x 1)     COMMANDED control [v; omega] (filter never sees the
%                       actuation noise -> genuine process uncertainty)
%   z       (2 x L)     measurements [range; bearing] from the TRUE pose
%   linit   (1 x L)     logical, whether each landmark is initialised
%   P       struct      .dt .L .M (2x2 control-noise cov) .Q (2x2 meas-noise cov)
%   do_update bool      true -> EKF-SLAM (predict+update); false -> dead reckoning
%
% The pure predict-only path (do_update=false) is the odometry baseline.

dt = P.dt; L = P.L;
n  = 3 + 2*L;
v  = u_cmd(1); w = u_cmd(2);
if abs(w) < 1e-9, w = 1e-9; end          % guard straight-line singularity
th = X(3);

F_x = [eye(3), zeros(3,2*L)];

% --- prediction (exact unicycle integration) ---
motion = [ v/w*( sin(th+w*dt) - sin(th) );
           v/w*( cos(th) - cos(th+w*dt) );
           w*dt ];
G_low = [0 0  v/w*( cos(th+w*dt) - cos(th) );
         0 0  v/w*( sin(th+w*dt) - sin(th) );
         0 0  0];
% control -> state noise mapping
V = [ ( -sin(th) + sin(th+w*dt) )/w,  v*( sin(th) - sin(th+w*dt) )/w^2 + v*cos(th+w*dt)*dt/w;
      (  cos(th) - cos(th+w*dt) )/w, -v*( cos(th) - cos(th+w*dt) )/w^2 + v*sin(th+w*dt)*dt/w;
        0,                             dt ];
R = V*P.M*V';

X = X + F_x'*motion;
G = eye(n) + F_x'*G_low*F_x;
Sigma = G*Sigma*G' + F_x'*R*F_x;

if ~do_update
    return            % odometry / dead-reckoning baseline: no correction
end

% --- sequential measurement update, one landmark at a time ---
for i = 1:L
    if ~linit(i)
        % initialise landmark from first observation (inverse measurement model)
        X(2+2*i) = X(1) + z(1,i)*cos(X(3)+z(2,i));
        X(3+2*i) = X(2) + z(1,i)*sin(X(3)+z(2,i));
        linit(i) = true;
    end
    dx = X(2+2*i) - X(1);
    dy = X(3+2*i) - X(2);
    q  = dx^2 + dy^2;
    sq = sqrt(q);
    z_hat = [sq; atan2(dy,dx) - X(3)];

    F_xj = zeros(5, n);
    F_xj(1:3,1:3) = eye(3);
    F_xj(4:5, 2+2*i:3+2*i) = eye(2);

    H = (1/q)*[ -sq*dx, -sq*dy,  0,  sq*dx,  sq*dy;
                  dy,    -dx,   -q,   -dy,    dx ] * F_xj;

    S = H*Sigma*H' + P.Q;
    K = Sigma*H'/S;

    innov = z(:,i) - z_hat;
    innov(2) = wrapToPi(innov(2));        % wrap bearing innovation
    X = X + K*innov;
    Sigma = (eye(n) - K*H)*Sigma;
end
end
