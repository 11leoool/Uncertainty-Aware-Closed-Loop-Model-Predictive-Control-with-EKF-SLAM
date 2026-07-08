function [o, So, oinit] = cv_tracker_um(o, So, z_o, x_est, Ppose, oinit, Pt)
% CV_TRACKER_UM  Constant-velocity EKF tracking a moving obstacle from
% range-bearing detections, unknown-map version:
%   * the obstacle is UNKNOWN until first detected (finite sensing range);
%   * initialisation inherits the robot pose uncertainty (like static_obs_ekf);
%   * pass z_o = [] for a step with NO detection (out of range): the tracker
%     then PREDICTS ONLY and its covariance grows -- which the covariance-aware
%     margin sees.
%
%   o   (4x1) [ox; oy; vx; vy],  So (4x4),  x_est robot pose estimate (3x1),
%   Ppose (3x3) its covariance,  Pt: .dt .sa2 .Ro

dt = Pt.dt;
F = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1];
G = [dt^2/2 0; 0 dt^2/2; dt 0; 0 dt];
Qo = G * (Pt.sa2*eye(2)) * G';

% --- initialise on first detection ---
if ~oinit
    if isempty(z_o)
        return                       % nothing to do until first detection
    end
    a = x_est(3) + z_o(2);
    o = [ x_est(1) + z_o(1)*cos(a);
          x_est(2) + z_o(1)*sin(a);
          0; 0 ];
    Gx = [1 0 -z_o(1)*sin(a);
          0 1  z_o(1)*cos(a)];
    Gz = [cos(a) -z_o(1)*sin(a);
          sin(a)  z_o(1)*cos(a)];
    So = blkdiag(Gx*Ppose*Gx' + Gz*Pt.Ro*Gz', 1.0*eye(2));  % unknown velocity
    oinit = true;
    return
end

% --- CV prediction ---
o  = F*o;
So = F*So*F' + Qo;

if isempty(z_o)
    return                           % out of range: prediction only
end

% --- measurement update about the estimated robot pose ---
dx = o(1) - x_est(1);
dy = o(2) - x_est(2);
q  = dx^2 + dy^2;  r = sqrt(q);
z_hat = [r; wrapToPi(atan2(dy,dx) - x_est(3))];
H = [ dx/r,  dy/r, 0, 0;
     -dy/q,  dx/q, 0, 0];
Hp = [ -dx/r, -dy/r,  0;
        dy/q, -dx/q, -1];            % d h / d robot pose
Reff = Pt.Ro + Hp*Ppose*Hp';         % fold pose uncertainty into meas noise
S = H*So*H' + Reff;
K = So*H'/S;
innov = z_o - z_hat;
innov(2) = wrapToPi(innov(2));
o  = o + K*innov;
So = (eye(4) - K*H)*So;
end
