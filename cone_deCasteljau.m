function [p,r_p] = cone_deCasteljau(t,x0,r_x0,v_x0,v_rx0,x3,r_x3,v_x3,v_rx3,delta)
% Conic de Casteljau from velocitiies at each point
% t ~ time between [0,1]
% x0,x3 ~ support knots
% r_x0, r_x3 ~ value knots
% v_x0,v_x3 ~ velocities of support knots
% v_rx0,v_rx3 ~ velocities of value at knots
% TODO take care of case s/r+3<0
n_v_x0 = sqrt(sum(v_x0.^2,2));
n_v_x3 = sqrt(sum(v_x3.^2,2));
c1 = max(atan( n_v_x0 ./ (v_rx0./r_x0 + 3/delta) ), 1e-3);

c2 = delta/3 * sqrt(n_v_x0.^2 + (v_rx0./r_x0+3/delta).^2);
c3 = max(atan(  n_v_x3 ./ (3/delta - v_rx3./r_x3) ), 1e-4);
c4 = delta/3 * sqrt(n_v_x3.^2 + (3/delta - v_rx3./r_x3).^2);

x1 = x0 + c1 .* v_x0 ./ n_v_x0;
r_x1 = r_x0 .* c2;
x2 = x3 - c3 .* v_x3 ./ n_v_x3;
r_x2 = r_x3 .* c4;

[w0,r_w0] = cone_geodesic2(t,x0,r_x0,x1,r_x1);
[w1,r_w1] = cone_geodesic2(t,x1,r_x1,x2,r_x2);
[w2,r_w2] = cone_geodesic2(t,x2,r_x2,x3,r_x3);

[u0,r_u0] = cone_geodesic2(t,w0,r_w0,w1,r_w1);
[u1,r_u1] = cone_geodesic2(t,w1,r_w1,w2,r_w2);

[p,r_p] = cone_geodesic2(t,u0,r_u0,u1,r_u1);
end