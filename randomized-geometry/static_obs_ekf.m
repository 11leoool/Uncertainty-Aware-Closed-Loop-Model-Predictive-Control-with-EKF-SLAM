function [o, Po] = static_obs_ekf(action, o, Po, z, xhat, Ppose, Rz)
% STATIC_OBS_EKF  2-state EKF for a STATIC obstacle sensed by range-bearing
% detections, consistent with the unknown-map setting: the obstacle is
% unknown until first detected, then estimated online from the robot's OWN
% (uncertain) pose estimate.
%
%   'init'   : initialise from the first detection (inverse measurement model);
%              Po = G_x P_pose G_x' + G_z R G_z'  -- the initial obstacle
%              covariance inherits the robot pose uncertainty.
%   'update' : measurement update. The robot pose enters the predicted
%              measurement; its uncertainty is folded into an effective
%              measurement covariance R_eff = R + H_pose P_pose H_pose'
%              (first-order), keeping Po honest without a joint filter.
%
%   z     : [range; bearing] detection
%   xhat  : robot pose estimate (3x1), Ppose its 3x3 covariance
%   Rz    : 2x2 detection noise covariance

switch action
    case 'init'
        a  = xhat(3) + z(2);
        o  = xhat(1:2) + z(1)*[cos(a); sin(a)];
        Gx = [1 0 -z(1)*sin(a);
              0 1  z(1)*cos(a)];
        Gz = [cos(a)  -z(1)*sin(a);
              sin(a)   z(1)*cos(a)];
        Po = Gx*Ppose*Gx' + Gz*Rz*Gz';

    case 'update'
        Po = Po + 1e-6*eye(2);                 % tiny process noise (static)
        dx = o(1) - xhat(1);  dy = o(2) - xhat(2);
        q  = dx^2 + dy^2;  sq = sqrt(q);
        z_hat = [sq; wrapToPi(atan2(dy,dx) - xhat(3))];
        Ho = (1/q)*[ sq*dx,  sq*dy;
                       -dy,     dx ];          % d h / d obstacle
        Hp = (1/q)*[ -sq*dx, -sq*dy,  0;
                        dy,    -dx,  -q ];     % d h / d robot pose
        Reff = Rz + Hp*Ppose*Hp';
        S = Ho*Po*Ho' + Reff;
        K = Po*Ho'/S;
        innov = z - z_hat;
        innov(2) = wrapToPi(innov(2));
        o  = o + K*innov;
        Po = (eye(2) - K*Ho)*Po;

    otherwise
        error('static_obs_ekf: unknown action %s', action);
end
end
