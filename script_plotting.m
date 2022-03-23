%% 1D Gaussians
clc, clear
p.n = 2^8; 
p.B = 0.06;
p.R = 2.7; % sqrt(-2*log(1e-2))

t = linspace(0,10,4); 
t = [0,10/3,2*10/3,10];
t = t(:);

mm = .5;
x{1} = linspace(mm-p.R*p.B,mm+p.R*p.B,p.n); x{1} = x{1}(:);
x{2} = [linspace(.3-p.R*p.B,.3+p.R*p.B,p.n) linspace(.7-p.R*p.B,.7+p.R*p.B,p.n)]; x{2} = x{2}(:);
x{3} = [linspace(.3-p.R*p.B,.3+p.R*p.B,p.n) linspace(.7-p.R*p.B,.7+p.R*p.B,p.n)]; x{3} = x{3}(:);
x{4} = linspace(0,1,p.n); x{4} = x{4}(:);

mu{1} = 1*exp(-.5*(x{1}-mm).^2/p.B^2);
mu{2} = .5*exp(-.5*(x{2}-.3).^2/p.B^2) + .5*exp(-.5*(x{2}-.7).^2/p.B^2);
mu{3} = exp(-.5*(x{3}-.3).^2/p.B^2) + exp(-.5*(x{3}-.7).^2/p.B^2);
mu{4} = 0*x{4} + .5;

tt = linspace(0,10,128); tt = tt(:);
B2 = .9;
xx = linspace(mm-B2*p.R*p.B,mm+B2*p.B*p.R,500); xx = xx(:);

p.UOT_eta = 50e2;
p.UOT_iter = 2e2;
p.plan2map_varianceQuantile = .9;
p.mins = 3;

[p_x,p_r] = WFR_spline(t,x,mu,tt,xx,p);

% Plotting
for ii = 1:1
fig1 = figure(1); clc
for i = 1:length(t)
    plot3(x{i},t(i)*ones(length(mu{i}),1),mu{i},'.-'), hold on
end

for i = 1:size(p_x,1)
    plot3(permute(p_x(i,:,:), [3 2 1]), tt, permute(p_r(i,:,:).^2, [3 2 1]), 'color',[.1,.1,.1,.15])
end
hold off
xlim([0 1]); ylim([0 10]); pbaspect([1 2 1]), grid on
zlim([0,2.5])
view(50,40)

xlabel('$x$'), ylabel('$t$')
legend({'\mu_1','\mu_2','\mu_3','\mu_4'}, ...
    'Location','northwest','NumColumns',2)
end


%%  2D Gaussians
clc, clear
p.n = 16;
p.B = .1;
p.UOT_eta = 10e3;
p.UOT_iter = 2e2;
p.plan2map_varianceQuantile = .8;
p.mins = 3;

t = linspace(0,1,4);
tt = linspace(0,1,64); tt = tt(:);

[xx,~] = gaussian_bump([0,0],5, 1.75*p.B);
xx = xx/10;

[x{1}, mu{1}] = gaussian_bump([0,0],p.n, 2*p.B);
mu{1} = .75*mu{1};

[dum1,dum2] = gaussian_bump(sqrt(2)/2 * [1,1], p.n, p.B);
[dum3,dum4] = gaussian_bump(sqrt(2)/2 * [1,-1],p.n,p.B);
[dum5,dum6] = gaussian_bump([1,0],p.n,p.B);
x{2} = [dum1;dum3;dum5];
mu{2} = .65*[dum2;dum4;dum6];

[dum1,dum2] = gaussian_bump([1.5,1], p.n, p.B);
[dum3,dum4] = gaussian_bump([1.5,-1],p.n,p.B);
x{3} = [dum1;dum3];
mu{3} = .75*[dum2;dum4];

[x{4}, mu{4}] = gaussian_bump([2,0],p.n, 2*p.B);

for i = 1:length(x); x{i} = x{i}/10; 
%     mu{i} = 0*mu{i} + 1;
end

[p_x,p_r] = WFR_spline(t,x,mu,tt,xx,p);

% variance of paths
clc
dx =  permute(p_x, [3 1 2]);
mux =  permute(p_r,[3 1 2]) ;
mux = mux ./ sum(mux,1);

meandx = sum(dx .* mux,1);

vardx = sum((dx - meandx).^2,[1 3]);
idx = find(vardx > quantile(vardx,.8));
 
% idx = unique([idx, [22 37 12 7 9 19 54 34] - 4]);


clear dum1 dum2 dum3 dum4 dum5 dum6

% Plotting fig1
for ii = 1:1
fig1 = figure(1); clc
for i = 1:4
%     ci = colmap(i,[0,1]);
%     ci = ci(2,:);
%     scatter3(x{i}(:,1),x{i}(:,2),1*mu{i},10,ci,"filled"), hold on 
    
% scatter3(x{i}(:,1),x{i}(:,2),mu{i}.^2,30,colmap(i,mu{i}.^2),"filled"), hold on  
    scatter3(x{i}(:,1),x{i}(:,2),mu{i}.^2,30,colmap(i,mu{i}.^2),"filled",'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2), hold on 

%     scatter(x{i}(:,1),x{i}(:,2),20,colmap(i,mu{i}),"filled"), hold on  % 2d
end
scatter3(p_x(:,1,1),xx(:,2,1),0*p_r(:,1,1),5,'red')
% scatter(p_x(:,1,1),xx(:,2,1),1,'red') % 2d

% %%
hold off, view(-40,65),


end

% Plotting fig2
for ii = 1:1

fig2 = figure(2); clc
for i = 1:4
    scatter(x{i}(:,1),x{i}(:,2),20,colmap(i,mu{i}),"filled"), hold on
%     scatter3(x{i}(:,1),x{i}(:,2),mu{i},30,colmap(i,mu{i}),"filled",'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1), hold on 
%     scatter3(x{i}(:,1),x{i}(:,2),mu{i},30,colmap(i,mu{i}),"filled"), hold on 
end

cd = colors(tt);
colormap(cd(1:3,:)')

for i = 1:size(p_x,1)
    pp = plot(permute(p_x(i,1,:),[3 1 2]), permute(p_x(i,2,:),[3 1 2]),'LineWidth',2 );
%     pp = plot3(permute(p_x(i,1,:),[3 1 2]), permute(p_x(i,2,:),[3 1 2]),permute(p_r(i,1,:),[3 1 2]).^2,'LineWidth',2 );
    drawnow
    set(pp.Edge, 'ColorBinding','interpolated', 'ColorData',cd);
end
scatter(p_x(:,1,1),xx(:,2,1),15,'red') % 2d

% %%
hold off, grid on
cb = colorbar; cb.Label.String = 'time';
xlabel("$x$"), ylabel('$y$')

end




% Plotting fig3
for ii = 1:1

fig3 = figure(3);
for i = 1:4
    scatter(x{i}(:,1),x{i}(:,2),20,colmap(i,mu{i}),"filled"), hold on
end

cd = colors(tt);

% for i = 1:max(idx)
maxcols = max(p_r(:))^2;
mincols = min(p_r(:))^2;
% mincols=0;
for i = idx
    cols2 = p_r(i,1,:).^2; cols2 = cols2(:);
    cols2 = 10*(cols2-mincols+.01)/maxcols;
    
    for ii = 1:(size(p_x,3)-1)
        plot(permute(p_x(i,1,ii + (0:1)),[3 1 2]), permute(p_x(i,2,ii + (0:1)),[3 1 2]), ...
            'LineWidth',cols2(ii),'color',cd(1:3,ii)');
        hold on
    end
end

scatter(p_x(idx,1,1),xx(idx,2,1),10,'red') % 2d

hold off, grid on

colormap(cd(1:3,:)')
cb = colorbar; hold off
cb.Label.String = 'time';
xlabel("$x$"), ylabel('$y$')


end



%% Aux functions for plotting
function cm = colmap(i,v)
% Helper function that returns color range 
% i ~ [1,2,3,4] choosing which color from the first four default
% v ~ values to be interpolated from max and 0
cc = [0, 0.4470, 0.7410;...
0.8500, 0.3250, 0.0980;...
0.9290, 0.6940, 0.1250;...
0.4940, 0.1840, 0.5560];
% cm = [flip(linspace(cc(1,i),1,n));
%     flip(linspace(cc(2,i),1,n));
%     flip(linspace(cc(3,i),1,n))];
% cm = cm';
v = v(:);
m = min(v);
m = 0; % TODO replace min by 0?
M = max(v);
cm = ( cc(i,:) .* (v-m) + [1,1,1] .* (M-v))/ (M-m);
end

function cols = colors(t)
% cols is 4 x length(t)
n = length(t);
cc = [0, 0.4470, 0.7410;...
0.8500, 0.3250, 0.0980;...
0.9290, 0.6940, 0.1250;...
0.4940, 0.1840, 0.5560]';
tt = linspace(0,1,4);
cols = zeros(4, n);
for i = 1:n
    idx = find(tt <= t(i), 1, 'last');
    if idx < 4
        cols(1:3,i) = (cc(:,idx) * (tt(idx+1)-t(i)) + cc(:,idx+1) * (t(i)-tt(idx))) / (tt(idx+1)-tt(idx));
    else
        cols(1:3,i) = (cc(:,idx-1) * (tt(idx)-t(i)) + cc(:,idx) * (t(i)-tt(idx-1))) / (tt(idx)-tt(idx-1));
    end
end

cols(1:3,:) = uint8(255*cols(1:3,:));
cols(4,:) = 1;
cols = uint8(cols);
end

function [x,fx] = gaussian_bump(center, n,B)
% Points sampled equispaced in a 2d sphere and their corresponding Gaussian
% density
% center ~ 2d vector with the mean
% n ~ number of radial points
% B ~ 
x = circle_supp(B*1.5,n) + center;
fx = exp(-.5*sum((x-center).^2,2)/B^2);
end


