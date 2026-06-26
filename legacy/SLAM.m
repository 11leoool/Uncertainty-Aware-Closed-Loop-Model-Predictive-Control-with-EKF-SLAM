function [x_output,y_output,theta_output,Sigma,X_p,estimate_location,Measurements] = SLAM(u_cl,t0,Measurements,Sigma,X_p,estimate_location,T)

delta_t = T;

kk = round(t0/delta_t+1);

covariance_omega = 0.000001;
covariance_v = 0.000001;
covariance_d = 0.01;
covariance_a = 0.01;
x_location = [-0.5 1 -1.5];
y_location = [0.5 1 0];

landmarks = length(x_location);


% states: x,y,theta t-1 state                   x(t-1) y(t-1) theta(t-1)
%controls: epsilon_v, epsilon_omega             u_cl(t-1,1);u_cl(t-1,2)
%current states:                                x(t) y(t) theta(t)
%measurements: distance, angle,                 Measurements
%covariance matrix: Sigma: size of (2N+3)*(2N+3)
%Noise in control space:                        covariance_omega,covariance_v
%Noise in measurements:                         covariance_a,covariance_d
%Noise in states:                               R_t
%Noise in measurements:                         Q_t



if kk == 1
    Sigma = [[zeros(3,3),zeros(3,2*landmarks)];[zeros(2*landmarks,3),1000000000000000*eye(2*landmarks,2*landmarks)]];
    X_p = [0;0;0;0;0;0;0;0;0];
end

x0 = [X_p(1);X_p(2);X_p(3)];

[Measurements] = measure(x0,t0,Measurements,T);

% Sigma = [zeros(9,9)+0.001]
R = [];
Q = [covariance_d 0;0 covariance_a];


F_x = [eye(3,3) zeros(3,2*landmarks)];




%
% X_p = [x;y;theta];

X_p_output(:,:,1) = X_p;
%for test
% load("matlab.mat")
% load("matlab_measurements_noisy_v_omega_d_a.mat")
% load("delta_x.mat")
% load("delta_y.mat")
Kalman_gain_total = [zeros(9,1)];
Kalman_gain_sigma_total = [zeros(9,9)];




epsilon_v = u_cl(round(kk-1),1);
epsilon_omega = u_cl(round(kk-1),2);

x(kk-1) = X_p(1);
y(kk-1) = X_p(2);
theta(kk-1) = wrapTo2Pi(X_p(3));

V_t1 = (-sin(theta(kk-1)) + sin(theta(kk-1)+epsilon_omega*delta_t))/epsilon_omega;
V_t2 = epsilon_v* (sin(theta(kk-1)) - sin(theta(kk-1)+epsilon_omega*delta_t))/(epsilon_omega)^2 ...
    + epsilon_v* cos(theta(kk-1)+epsilon_omega*delta_t)*delta_t/epsilon_omega;
V_t3 = (cos(theta(kk-1)) - cos(theta(kk-1)+epsilon_omega*delta_t))/epsilon_omega;
V_t4 = -epsilon_v* (cos(theta(kk-1)) - cos(theta(kk-1)+epsilon_omega*delta_t))/(epsilon_omega)^2 ...
    + epsilon_v* sin(theta(kk-1)+epsilon_omega*delta_t)*delta_t/epsilon_omega;
V_t5 = 0;
V_t6 = delta_t;

M_t = [covariance_v 0;0 covariance_omega];

V_t = [V_t1 V_t2;
    V_t3 V_t4;
    V_t5 V_t6];

R = V_t*M_t*V_t';



%     theta(t) = wrapTo2Pi(theta(t-1) + epsilon_omega*delta_t);


motion_model = [-epsilon_v/epsilon_omega*sin(theta(kk-1)) + epsilon_v/epsilon_omega* sin(theta(kk-1) + epsilon_omega*delta_t);
    epsilon_v/epsilon_omega*cos(theta(kk-1)) - epsilon_v/epsilon_omega* cos(theta(kk-1) + epsilon_omega*delta_t);
    epsilon_omega*delta_t];


motion_model_Jacobian = [0 0 -epsilon_v/epsilon_omega*cos(theta(kk-1)) + epsilon_v/epsilon_omega* cos(theta(kk-1) + epsilon_omega*delta_t);
    0 0 -epsilon_v/epsilon_omega*sin(theta(kk-1)) + epsilon_v/epsilon_omega* sin(theta(kk-1) + epsilon_omega*delta_t);
    0 0 0];

X_p = X_p + F_x' * motion_model;
G_t = eye(9,9) + F_x' * motion_model_Jacobian * F_x;
Sigma = G_t*Sigma*G_t' + F_x'*R*F_x;


F_xj(:,:,1) = [eye(3,3) zeros(3,6);zeros(2,3) eye(2,2) zeros(2,4)];
F_xj(:,:,2) = [eye(3,3) zeros(3,6);zeros(2,3) zeros(2,2) eye(2,2) zeros(2,2)];
F_xj(:,:,3) = [eye(3,3) zeros(3,6);zeros(2,3) zeros(2,4) eye(2,2)];
for i = 1:landmarks
    if kk == 2
        % First observation: initialise the landmark estimate IN THE STATE
        % vector X_p (slots 2+2i, 3+2i) from the inverse measurement model.
        X_p(2+2*i) = X_p(1) + Measurements(2*kk-1,i)*cos(X_p(3)+Measurements(2*kk,i));
        X_p(3+2*i) = X_p(2) + Measurements(2*kk-1,i)*sin(X_p(3)+Measurements(2*kk,i));
    end
    % Predicted measurement is computed from the ESTIMATED landmark state
    % (true SLAM: no ground-truth landmark coordinates are used here).
    delta_x(i) = X_p(2+2*i) - X_p(1);
    delta_y(i) = X_p(3+2*i) - X_p(2);
    q(i) = delta_x(i)^2 + delta_y(i)^2;
    z_hat_d(i) = sqrt(q(i));
    z_hat_a(i) = atan2(delta_y(i), delta_x(i)) - X_p(3);
    % keep the returned helper variable consistent with the state estimate
    estimate_location(:,:,i) = [X_p(2+2*i); X_p(3+2*i)];
end

for i = 1:landmarks
    H(:,:,i) = 1/q(i)*[-sqrt(q(i))*delta_x(i), -sqrt(q(i))*delta_y(i), 0, sqrt(q(i))*delta_x(i), sqrt(q(i))*delta_y(i);
        delta_y(i), -delta_x(i), -q(i), -delta_y(i), delta_x(i)]*F_xj(:,:,i);
end


for i = 1:landmarks
    K(:,:,i) = Sigma*H(:,:,i)'* inv(H(:,:,i)*Sigma*H(:,:,i)'+Q);

end

for i = 1:landmarks
    innov = [Measurements(2*kk-1,i);Measurements(2*kk,i)] - [z_hat_d(i);z_hat_a(i)];
    innov(2) = wrapToPi(innov(2));   % wrap bearing innovation to [-pi,pi]
    Kalman_gain(:,:,i) = K(:,:,i)*innov;
    Kalman_gain_total = Kalman_gain_total + Kalman_gain(:,:,i);


    Kalman_gain_sigma(:,:,i) = K(:,:,i)*H(:,:,i);
    Kalman_gain_sigma_total = Kalman_gain_sigma_total + Kalman_gain_sigma(:,:,i);

end

X_p = Kalman_gain_total + X_p;


Sigma = (eye(9,9)-Kalman_gain_sigma_total)*Sigma;



Kalman_gain_total = [zeros(9,1)];
Kalman_gain_sigma_total = [zeros(9,9)];

X_p;


x_output = X_p(1);
y_output = X_p(2);
theta_output = X_p(3);
%     X_p_output(:,:,t) = X_p;



clc



end





