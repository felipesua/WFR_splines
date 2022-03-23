function [p,r_p] = WFR_spline(t,x,mu,tt,xx,p)
% Wasserstein-Fisher-Rao spline 
% x ~ (n_pts, dims, n_times) tensor of supports of knots (measures)
% mu ~ (n_pts, 1, n_times) tensor of values of measures
% tt ~ interpolating times 
% xx ~ initial (at time 0) interpolating support 
% p parameters of algorithm:
%     p.UOT_eta  Unbalanced OT regularization 
%     p.UOT_iter Unbalanced OT iterations
%     p.plan2map_varianceQuantile quantile of variance for plan to map projection
%     p.mins number of min distance to construct interpolation in support

assert(issorted(t),'time of knots is not sorted')
assert(issorted(tt), 'interpolation time is not sorted')
assert(length(x) == length(mu), 'x and mu(x) have different lengths')
assert(length(t) == length(x), 'x and t must have different lengths')
assert(all(diff(cellfun(@(a) size(a,2),x)) == 0),'not all x have the same dimension')
assert((t(1) <= tt(1)) & (t(end) >= tt(end)), 'interpolating times lie outside of knot times.')

if ~isfield(p,'UOT_eta')
    p.UOT_eta = 5e3;
end
if ~isfield(p,'UOT_iter')
    p.UOT_iter = 2e2;
end
if ~isfield(p,'plan2map_varianceQuantile')
    p.plan2map_varianceQuantile = .8;
end

n_knots = length(x);
d = size(x{1},2);

% Plans M
M = cell(n_knots-1,1);

for i = 1:(n_knots-1)
    C = pdist2(x{i},x{i+1});
    assert(max(C(:)) <=  3.141593/2,'Distances exceed pi/2');
    M{i} = UOT(p.UOT_eta, -2*log(cos(C)),  mu{i}, mu{i+1}, p.UOT_iter);
end

% map knots (xi,yi)
xi = cell(n_knots-1,1); fxi = cell(n_knots-1,1);
yi = cell(n_knots-1,1); fyi = cell(n_knots-1,1);
for i = 1:(n_knots-1)
    [xi{i}, fxi{i}, yi{i}, fyi{i}] = plan2map(M{i},x{i},mu{i},x{i+1},mu{i+1},p.plan2map_varianceQuantile);
end

% mapped points (xi_,fxi)_j --> (xi_,fxi_)_{j+1}
n_samples = size(xx,1);
xi_ = zeros([n_samples, d, n_knots]);
fxi_ = zeros([n_samples, 1, n_knots]);

if ~isfield(p,'mins')
    p.mins = 1.5;
end

xi_(:,:,1) = xx;

for i = 1:(n_knots-1)
    W = min_interp(xi{i},xi_(:,:,i),p.mins,d);
    fxi_(:,1,i) = W * fxi{i} ./ sum(W,2);
    xi_(:,:,i+1) = W * yi{i} ./ sum(W,2);
end

fxi_(:,1,n_knots) = W * fyi{n_knots-1} ./ sum(W,2);

% derivatives
vx = tensor_spline_derivative(t,xi_);
vr = tensor_spline_derivative(t,sqrt(fxi_));

% Interpolation
n_t = length(tt);

p = zeros(n_samples,d,n_t);
r_p = zeros(n_samples,1,n_t);

for i = 1:n_t
    if t(end) <= tt(i)
        idx = n_knots-1;
    else
        idx = find(t > tt(i),1) - 1;
    end
    delta = t(idx+1) - t(idx);
    [pi,r_pi] = cone_deCasteljau( (tt(i) - t(idx))/delta, ...
        xi_(:,:,idx), sqrt(fxi_(:,:,idx)),  vx(:,:,idx), vr(:,:,idx),...
        xi_(:,:,idx+1), sqrt(fxi_(:,:,idx+1)),  vx(:,:,idx+1), vr(:,:,idx+1),...
        delta);
    p(:,:,i) = pi;
    r_p(:,:,i) = r_pi;
end

end