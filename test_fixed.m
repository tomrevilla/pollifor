% This script generates time series for fixed preferences
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
% pref = 0: REPLICATOR DYNAMICS. Here u is drawn from the [0,1] interval and 
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
c = 0.4*[1, 1];
K = 50*[1, 1];
g = 0.0;     

% Operational parameters:
options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6 1e-6], 'MaxStep', 1e-10);
tspan = 100000;
pref = 0;


% LEFT-UP: Invasion and equilibrium
K = 60*[1, 1];
c = 0.4*[1, 1];
u = 0.8;
y0 = [43.5210, 0.01, 31.6893];
tspan = 2000;
[t, y] = ode45(@odesystem, [0, tspan], y0, [], r, e, [u, 1-u].*b, a, w, m, c, K, g, pref);
y = max(0,y);
yf1 = y(length(y)-round(length(y)/4):end,:);

subplot(2,2,1)
hold on
plot(t,y(:,1), 'g','LineWidth',1)
plot(t,y(:,2), 'r','LineWidth',1)
plot(t,y(:,3), 'm','LineWidth',1)
plot(t(1:20:end),y(1:20:end,1), 'gs','MarkerSize',4)
plot(t(1:20:end),y(1:20:end,2), 'rd','MarkerSize',4)
plot(t(1:20:end),y(1:20:end,3), 'mo','MarkerSize',4)
axis([0 tspan 0 50])
set(gca,'YTick',[10:10:50])
box on
%xlabel('time','fontsize',12)
ylabel('P_1 , P_2 , A','fontsize',12)
title('(a)','fontsize',12)
clear y u

% RIGHT-UP: Limit cycles
K = 60*[1, 1];
c = 1.2*[1, 1];
u = 0.605;
y0 = [50, 100, 50];
tspan = 40000;
[t, y] = ode45(@odesystem, [0, tspan], y0, [], r, e, [u, 1-u].*b, a, w, m, c, K, g, pref);
y = max(0,y);
yf3 = y(length(y)-round(length(y)/4):end,:);

subplot(2,2,2)
hold on
plot(t,y(:,1), 'g','LineWidth',1)
plot(t,y(:,2), 'r','LineWidth',1)
plot(t,y(:,3), 'm','LineWidth',1)
plot(t(1:200:end),y(1:200:end,1), 'gs','MarkerSize',4)
plot(t(1:200:end),y(1:200:end,2), 'rd','MarkerSize',4)
plot(t(1:200:end),y(1:200:end,3), 'mo','MarkerSize',4)
axis([0 tspan 0 40])
set(gca,'YTick',[10:10:40])
box on
%xlabel('time','fontsize',12)
%ylabel('P_1 , P_2 , A','fontsize',12)
title('(b)','fontsize',12)
clear y u
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:))));


% LEFT-DOWN: Damped oscillations
K = 60*[1, 1];
c = 1.2*[1, 1];
u = 0.607;
y0 = [50, 100, 50];
tspan = 40000;
[t, y] = ode45(@odesystem, [0, tspan], y0, [], r, e, [u, 1-u].*b, a, w, m, c, K, g, pref);
y = max(0,y);
yf4 = y(length(y)-round(length(y)/4):end,:);

subplot(2,2,3)
hold on
plot(t,y(:,1), 'g','LineWidth',1)
plot(t,y(:,2), 'r','LineWidth',1)
plot(t,y(:,3), 'm','LineWidth',1)
plot(t(1:200:end),y(1:200:end,1), 'gs','MarkerSize',4)
plot(t(1:200:end),y(1:200:end,2), 'rd','MarkerSize',4)
plot(t(1:200:end),y(1:200:end,3), 'mo','MarkerSize',4)
axis([0 tspan 0 40])
set(gca,'YTick',[10:10:40])
box on
xlabel('time','fontsize',12)
ylabel('P_1 , P_2 , A','fontsize',12)
title('(c)','fontsize',12)
clear y u
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:))));


% RIGHT-DOWN: Invasion and collapse
K = 35*[1, 1];
c = 0.8*[1, 1];
u = 0.4;
y0 = [21.1086, 0.01, 10.6380];
tspan = 10000;
[t, y] = ode45(@odesystem, [0, tspan], y0, [], r, e, [u, 1-u].*b, a, w, m, c, K, g, pref);
y = max(0,y);
yf2 = y(length(y)-round(length(y)/4):end,:);

subplot(2,2,4)
hold on
plot(t,y(:,1), 'g','LineWidth',1)
plot(t,y(:,2), 'r','LineWidth',1)
plot(t,y(:,3), 'm','LineWidth',1)
plot(t(1:50:end),y(1:50:end,1), 'gs','MarkerSize',4)
plot(t(1:50:end),y(1:50:end,2), 'rd','MarkerSize',4)
plot(t(1:50:end),y(1:50:end,3), 'mo','MarkerSize',4)
axis([0 tspan 0 25])
set(gca,'YTick',[5:5:25])
box on
xlabel('time','fontsize',12)
%ylabel('P_1 , P_2 , A','fontsize',12)
title('(d)','fontsize',12)
clear y u


elapsed_time = toc

print -depsc test_fixed.eps
print  -dpng test_fixed.png

