function W = min_interp(x,x_,B,d)
% Interpolation using k-nn and at most d points
% x interpolating poitns
% x_ interpolated points
% B number of min distance to look for
% d dimension
W = pdist2(x_,x);
min_d = mink(W,d+2,2);

W = sparse(W <= min(B*min_d(:,1), min_d(:,end)));
end