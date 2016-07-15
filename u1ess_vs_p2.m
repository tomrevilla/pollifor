% Preference for focal plant against competitor
clf
clear

% Default parameters
e = 0.1*[2, 1];
b = 0.1*[1, 1];
a = 0.4*[1, 1];
w = 0.25*[1, 1];

P1 = 10;
P2 = linspace(0,50,200);
A = [0 2];

u1g = e(1)*a(1)*P1./(e(1)*a(1)*P1 + e(2)*a(2)*P2);

u1s = (w(2)*e(1)*a(1)*b(1)*P1 - w(1)*e(2)*a(2)*b(2)*P2)./(b(1)*b(2)*A(1)*(e(1)*a(1)*P1 + e(2)*a(2)*P2));
plot(P2, max(0, min(1, u1g + u1s)), '--k')
hold on

u1s = (w(2)*e(1)*a(1)*b(1)*P1 - w(1)*e(2)*a(2)*b(2)*P2)./(b(1)*b(2)*A(2)*(e(1)*a(1)*P1 + e(2)*a(2)*P2));
plot(P2, max(0, min(1, u1g + u1s)), 'k', 'linewidth', 2)

plot(P2, max(0, min(1, u1g)), 'k')

legend('A\rightarrow 0','A = 2','A\rightarrow \infty')
legend('boxoff')
axis([0 inf -0.1 1.1])
xlabel('Plant 2 density (P_2)','fontsize',12)
ylabel('ESS preference for plant 1 (u_1^*)','fontsize',12)

set(gca,'fontsize',12)
print -depsc u1ess_vs_p2.eps
print  -dpng u1ess_vs_p2.png