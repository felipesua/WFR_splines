function X = circle_supp(R,n)
% Create quasi-equispaced points in a 2d-ball 
% R ~ Radius of ball
% n ~ number of radial points
B = 4*n/R;
Rs = linspace(0,R,n);
x = [0]; y = [0];
for i = 1:n
x = [x Rs(i)*cos( linspace(0,2*pi,floor( B*Rs(i)) ) )];
y = [y Rs(i)*sin( linspace(0,2*pi,floor( B*Rs(i)) ) )];
end
X = [x(:) y(:)];
end