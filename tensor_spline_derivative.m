function v = tensor_spline_derivative(t,x)
% Tensor version of computation of derivatives at knots in spline
% t ~ knot_times
% x ~interpolating values (support or measure value)
t = t(:);
delta = diff(t);

T = 2*diag(delta(1:(end-1)) + delta(2:end)) + ...
    diag(delta(2:(end-1)),1) + diag(delta(2:(end-1)),-1);
[pts,dims,times] = size(x);

x = reshape(x, [pts*dims, times]);

Del = (x(:,3:end) - x(:,2:(end-1))) ./ delta(2:end)' - ...
    (x(:,2:(end-1)) - x(:,1:(end-2))) ./ delta(1:(end-1))' ;
m = 6 * Del * inv(T); % TODO test backslash instead of inv
m = [0*m(:,1) m 0*m(:,1)];
v = (x(:,2:end) - x(:,1:(end-1))) ./ delta' - 1/6 * (m(:,2:end) + 2*m(:,1:(end-1))) .* delta';

a = (m(:,end) - m(:,end-1)) / (6 * delta(end));
b = m(:,end-1)/2;

vf = 3*delta(end)^2 * a + 2*delta(end) * b + v(:,end);
v = [v vf];

v = reshape(v,[pts,dims,times]);
end