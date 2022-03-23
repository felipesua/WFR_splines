
function [y,R] = cone_geodesic2(t,x0,r0,x1,r1)
% Cone geodesic at times t
% x0,x1 ~x endpoints
% r0, r1 ~r endpoints

t = permute(t(:), [2 3 1]);
th = sqrt(sum((x0-x1).^2,2));
cth = cos(th);
R = sqrt( (1-t).^2 .* (r0.^2) +t.^2 .* (r1.^2) + 2*(t.*(1-t)) .* (r0.*r1.*cth));
rho = acos( ((1-t) .* r0 + t .* (r1.*cth)) ./ R ) ./ th ;

rho(th < 1e-8,1,:) = t .* ones(sum(th < 1e-8),1);

y = x0 .* (1-rho) + x1 .* rho;
end

function X = circle_supp(R,n)
B = 4*n/R;
Rs = linspace(0,R,n);
x = [0]; y = [0];
for i = 1:n
x = [x Rs(i)*cos( linspace(0,2*pi,floor( B*Rs(i)) ) )];
y = [y Rs(i)*sin( linspace(0,2*pi,floor( B*Rs(i)) ) )];
end
X = [x(:) y(:)];
end