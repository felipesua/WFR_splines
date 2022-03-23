function [x_,fx_,y_,fy_] = plan2map(M,x,fx,y,fy,qt)
% Find a map from a plan
% M ~ plan
% x,y ~support of measures
% fx, fy ~ value of measures at support x,y
% qt ~ quantiles of variance to subset from
varis = zeros(1,size(M,1));

for i = 1:size(M,1)
s = sum(M(i,:));
varis(i) = sum((M(i,:) * y.^2 / s)  - (M(i,:) *y / s).^2);
end

varis2 = zeros(1,size(M,2));
for j = 1:size(M,2)
s = sum(M(:,j));
varis2(j) = sum((M(:,j)' * x.^2 / s)  - (M(:,j)' * x / s).^2);
end

q = quantile([varis, varis2], qt);
% size( (M(:,varis2<q)'*x) ./ sum(M(:,varis2<q),1)' )
x_ = [x(varis<q,:) ; (M(:,varis2<q)'*x) ./ sum(M(:,varis2<q),1)'];
y_ = [M(varis<q,:)*y ./ sum(M(varis<q,:),2); y(varis2<q,:)];
fx_ = [fx(varis<q) ; (M(:,varis2<q)'*fx) ./ sum(M(:,varis2<q),1)'];
fy_ = [M(varis<q,:)*fy ./ sum(M(varis<q,:),2); fy(varis2<q)];
end