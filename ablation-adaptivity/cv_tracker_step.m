function [o, So, oinit] = cv_tracker_step(o, So, z_o, x_est, oinit, Pt)
% CV_TRACKER_STEP  One predict+update of a constant-velocity (CV) EKF that tracks
% a moving obstacle from the robot's range-bearing detections.
%
%   o      (4x1)  obstacle state estimate [ox; oy; vx; vy]
%   So     (4x4)  covariance
%   z_o    (2x1)  range-bearing measurement of the obstacle [range; bearing]
%                 (generated from the TRUE robot pose; here interpreted via x_est)
%   x_est  (3x1)  robot pose ESTIMATE (the tracker only knows the SLAM estimate)
%   oinit  bool   whether the track is initialised
%   Pt     struct .dt, .sa2 (accel-noise variance), .Ro (2x2 meas-noise cov)
%
% Returns the posterior obstacle state/covariance. The position block So(1:2,1:2)
% (and its horizon roll-out) is the obstacle-prediction uncertainty used by the
% covariance-aware margin.

dt = Pt.dt;

% --- initialise on first detection (from estimated pose + relative measurement) ---
if ~oinit
    o = [ x_est(1) + z_o(1)*cos(x_est(3)+z_o(2));
          x_est(2) + z_o(1)*sin(x_est(3)+z_o(2));
          0; 0 ];
    So = diag([0.25, 0.25, 1.0, 1.0]);   % unknown initial velocity -> large vel cov
    oinit = true;
end

% --- CV prediction ---
F = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1];
G = [dt^2/2 0; 0 dt^2/2; dt 0; 0 dt];
Qo = G * (Pt.sa2*eye(2)) * G';
o  = F*o;
So = F*So*F' + Qo;

% --- range-bearing measurement update about the estimated robot pose ---
dx = o(1) - x_est(1);
dy = o(2) - x_est(2);
q  = dx^2 + dy^2;  r = sqrt(q);
z_hat = [r; wrapToPi(atan2(dy,dx) - x_est(3))];
H = [ dx/r,  dy/r, 0, 0;
     -dy/q,  dx/q, 0, 0];
S = H*So*H' + Pt.Ro;
K = So*H'/S;
innov = z_o - z_hat;
innov(2) = wrapToPi(innov(2));
o  = o + K*innov;
So = (eye(4) - K*H)*So;
end
