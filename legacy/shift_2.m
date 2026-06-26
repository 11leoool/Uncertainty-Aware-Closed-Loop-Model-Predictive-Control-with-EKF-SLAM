function [t0, x0, u0] = shift_2(T, t0, x0, u,f)
%with noise
st = x0;
con = u(1,:)';
% f_value = f(st,con)+ normrnd(0,0.01);

f_value = f(st,con);

st = st+ (T*f_value);
% st = st+ (T*f_value)+normrnd(0,0.0001);


x0 = full(st);

% t0 = t0 + T;
u0 = [u(2:size(u,1),:);u(size(u,1),:)];
end