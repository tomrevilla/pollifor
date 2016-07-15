% This script generates time series for adaptive preferences
%
% Variables and Parameters:
%
% y: vector of population densities and preference (*)
% r: vector of plant yields
% e: vector of animal yields
% b: vector of consumption (pollination or frugivory) rates (*)
% a: vector of (flower, fruit) supply rates
% w: vector of loss rates
% m: vector of plant and animal mortality rates
% c: vector of plant competition coefficients
% K: vector of plant zero yield density
% v: fixed or initial relative preference for plant 1 (1-v for plant 2) (*)
% g: adaptation rate in flexible diet models
%
% (*) Preference models
%
% pref = 0: FIXED PREFERENCES. Consumption rates are multiplied by v and 1-v and 
%           they remain so.
% pref = 1: ESS. Consumption rates are multiplied by u and 1-u, where u is the 
%           best response B(u) under current plant and pollinator densities.
% pref = 2: BEST RESPONSE DYNAMICS. Here u is drawn from the [0,1] interval and 
%           then follows the du/dt = g*(B(u) - u) differential equation
% pref = 3: REPLICATOR DYNAMICS. Here u is drawn from the [0,1] interval and 
%           then follows the replicator equation.
%

clear
clf
tic

% Default parameters, mostly symmetric
r = 0.1*[1, 1];
e = 0.1*[2, 1];
b = 0.1*[1, 1];
a = 0.4*[1, 1];
w = 0.25*[1, 1];
m = [0.01, 0.01, 0.1];
c = 0.0*[1, 1];
K = 50*[1, 1];

% Operational parameters:
options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6 1e-6], 'MaxStep', 1e-10);
pref = 3;

v = 0.5;
y0 = [v*K(1), (1-v)*K(2), 1, 0.999]
tspan = 10000;
pts=floor(logspace(0,log10(tspan)-1,15));
ptx = logspace(0,log10(tspan),log10(tspan)+1);

% Infinite
g = 0;
[t, y] = ode45(@odesystem, [0, tspan], y0(1:3), [], r, e, b, a, w, m, c, K, g, 1);
y = max(0,y);

for i = 1:length(t)
  u(i) = best_resp_u(y(i,1:3), e.*a, b, w);
end

subplot(2,2,1)
[ax, h1, h2]=plotyy(t,y(:,1:3),t,-1*ones(1,length(t)))
set(ax(1),'YLim',1.1*[0 50],'YTick',[10:10:50])
set(ax(2),'YLim',1.1*[0 1],'YTick',[0.2:0.2:1],'YColor','black')

hold(ax(1))
plot(t,y(:,1), 'g','LineWidth',1)
plot(t,y(:,2), 'r','LineWidth',1)
plot(t,y(:,3), 'm','LineWidth',1)

plot(t(pts),y(pts,1), 'gs','MarkerSize',4)
plot(t(pts),y(pts,2), 'rd','MarkerSize',4)
plot(t(pts),y(pts,3), 'mo','MarkerSize',4)

hold(ax(2))                      
%plot(ax(2),t,y(:,4), '--k','LineWidth',1)
plot(ax(2),t,u, '-k','LineWidth',.5)
plot(ax(2),t(pts),u(pts), '*k','MarkerSize',5)

title('(a: \nu=\infty)','fontsize',12)
%xlabel('time (log scale)','fontsize',12)
ylabel(ax(1),'P_1 , P_2 , A','fontsize',12) % right y-axis
%ylabel(ax(2),'   u_1','Rotation',0,'fontsize',12) % right y-axis

set(ax,'XScale','log');
set(ax,'Xtick',ptx);

clear y u g

% Fast
g = 1;
[t, y] = ode45(@odesystem, [0, tspan], y0, [], r, e, b, a, w, m, c, K, g, pref);
y = max(0,y);
yf3 = y(end,:)
for i = 1:length(t)
  u(i) = best_resp_u(y(i,1:3), e.*a, b, w);
end

subplot(2,2,2)
[ax, h1, h2]=plotyy(t,y(:,1:3),t,-1*ones(1,length(t)))
set(ax(1),'YLim',1.1*[0 50],'YTick',[10:10:50])
set(ax(2),'YLim',1.1*[0 1],'YTick',[0.2:0.2:1],'YColor','black')

hold(ax(1))
plot(t,y(:,1), 'g','LineWidth',1)
plot(t,y(:,2), 'r','LineWidth',1)
plot(t,y(:,3), 'm','LineWidth',1)

plot(t(pts),y(pts,1), 'gs','MarkerSize',4)
plot(t(pts),y(pts,2), 'rd','MarkerSize',4)
plot(t(pts),y(pts,3), 'mo','MarkerSize',4)

hold(ax(2))                      
plot(ax(2),t,y(:,4), 'k','LineWidth',.5)
%plot(ax(2),t,u, '-k','LineWidth',.5)
%plot(ax(2),t(pts),u(pts), '*k','MarkerSize',6)

title('(b: \nu=1)','fontsize',12)
xlabel('time (log scale)','fontsize',12)
%ylabel(ax(1),'P_1 , P_2 , A','fontsize',12) % right y-axis
ylabel(ax(2),'   u_1','Rotation',0,'fontsize',12) % right y-axis

set(ax,'XScale','log');
set(ax,'Xtick',ptx);

clear y u g

% Slow
g = 0.25;
[t, y] = ode45(@odesystem, [0, tspan], y0, [], r, e, b, a, w, m, c, K, g, pref);
y = max(0,y);

for i = 1:length(t)
  u(i) = best_resp_u(y(i,1:3), e.*a, b, w);
end

subplot(2,2,3)
[ax, h1, h2]=plotyy(t,y(:,1:3),t,-1*ones(1,length(t)))
set(ax(1),'YLim',1.1*[0 50],'YTick',[10:10:50])
set(ax(2),'YLim',1.1*[0 1],'YTick',[0.2:0.2:1],'YColor','black')

hold(ax(1))
plot(t,y(:,1), 'g','LineWidth',1)
plot(t,y(:,2), 'r','LineWidth',1)
plot(t,y(:,3), 'm','LineWidth',1)

plot(t(pts),y(pts,1), 'gs','MarkerSize',4)
plot(t(pts),y(pts,2), 'rd','MarkerSize',4)
plot(t(pts),y(pts,3), 'mo','MarkerSize',4)

hold(ax(2))                      
plot(ax(2),t,y(:,4), 'k','LineWidth',.5)
%plot(ax(2),t,u, '-k','LineWidth',.5)
%plot(ax(2),t(pts),u(pts), '*k','MarkerSize',6)

title('(c: \nu=0.25)','fontsize',12)
xlabel('time (log scale)','fontsize',12)
ylabel(ax(1),'P_1 , P_2 , A','fontsize',12) % right y-axis
ylabel(ax(2),'   u_1','Rotation',0,'fontsize',12) % right y-axis

set(ax,'XScale','log');
set(ax,'Xtick',ptx);

clear y u g


elapsed_time = toc

print -depsc test_replicator.eps
print  -dpng test_replicator.png

