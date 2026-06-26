function [Measurements,epsilon_a,epsilon_d] = measure(x0,t0,Measurements,T)
x = x0(1);
y = x0(2);
theta = x0(3);
delta_t = T;
kk = round(t0/delta_t+1);

x_location = [-0.5 1 -1.5];
y_location = [0.5 1 0];

% Measurement-noise variances used to GENERATE the (noisy) measurements.
% These must match the filter's assumed Q in SLAM.m (diag([0.01,0.01]))
% so the EKF is neither over- nor under-confident.
covariance_d = 0.01;   % range  noise variance  -> std 0.1 m
covariance_a = 0.01;   % bearing noise variance  -> std 0.1 rad (~5.7 deg)


landmarks = length(x_location);

%measurement
for i = 1:landmarks
    epsilon_d = mvnrnd(0,covariance_d);
    epsilon_a = mvnrnd(0,covariance_a);
    d(kk,i) = sqrt((x-x_location(i))^2+(y-y_location(i))^2) + epsilon_d;
    a(kk,i) = wrapTo2Pi(wrapTo2Pi(wrapTo2Pi(atan2((y_location(i)-y),(x_location(i)-x)))-wrapTo2Pi(theta))+epsilon_a);
end

M_store = [d(kk,:);a(kk,:)];
Measurements = [Measurements;M_store];




end