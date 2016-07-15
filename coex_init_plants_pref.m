% GRAPHICAL RESULTS -- SCENARIO 2
clf
clear

% Color codes: white ; green ; red   ; magenta
cmap =       [ 1,1,1 ; 0,1,0 ; 1,0,0 ; 1,0.5,1 ];
%cmap =       [ 1,1,1 ; 0.6,0.6,0.6 ; 0.8,0.8,0.8 ; 0.7,0.7,0.7 ];
colormap(cmap)

% Global axes properties
set(0,'defaultaxesfontsize',12);
equil_str = {'P_1 = ';'P_2 = ';'A = ';'u_1 = '};
alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p'];
tiKs = 0:0.2:1;

% Data strings. This is just an example. Adjust according to filename structure
row = 4;
col = 4;
datastring=[
'c_0.0/data_g0000.mat';
'c_0.4/data_g0000.mat';
'c_0.8/data_g0000.mat';
'c_1.2/data_g0000.mat';
'c_0.0/data_g0010.mat';
'c_0.4/data_g0010.mat';
'c_0.8/data_g0010.mat';
'c_1.2/data_g0010.mat';
'c_0.0/data_g0100.mat';
'c_0.4/data_g0100.mat';
'c_0.8/data_g0100.mat';
'c_1.2/data_g0100.mat';
'c_0.0/data_g9999.mat';
'c_0.4/data_g9999.mat';
'c_0.8/data_g9999.mat';
'c_1.2/data_g9999.mat'
];

% Get ESS preferences and additional settings from any file
load(datastring(1,:), 'gridres','init_pla','init_pol','K','e','a','b','w')
u = zeros(gridres);
maxpop = mean(K);
for xx = 1:gridres
  x(1:2) = [init_pla(xx), maxpop - init_pla(xx)];
  for yy = 1:gridres
    x(3) = init_pol;
    u(xx,yy) = best_resp_u(x, e.*a, b, w);
  end
end
tiXs = maxpop*tiKs;
tiYs = tiKs;
clear gridres init_pla init_pol maxpop e a b w

% For all plots
for panel=1:row*col
  % Load data and label outcomes
  load(datastring(panel,:), 'X','Y','Pop1','Pop2','Pop3','Pref','extinc_level','g','c')
  
  % Numeric codes:
  P1=1*(Pop1>extinc_level); % 0: P1 extinct, 1: P1 present
  P2=2*(Pop2>extinc_level); % 0: P2 extinct, 2: P2 present
  
  subplot(row,col,panel)

  % Coexistence and winner zones
  % if P1+P2 = 0: extinction
  % if P1+P2 = 1: species 1 wins
  % if P1+P2 = 2: species 2 wins
  % if P1+P2 = 3: coexistence
  imagesc(tiXs,tiYs,(P1+P2)')
  hold on
  
  % ESS 0 and 1 contour lines
%  contour(X,Y,u',[0,0],'k');
%  contour(X,Y,u',[1,1],'k');
  set(gca,'ticklength',2*get(gca,'ticklength'))
  % Plot labels
  if mod(panel-1,col) == 0
    set(gca,'YTick', tiYs,'YDir','normal', 'FontSize',10)
    ylabel('u_1(0)      ')
    set(get(gca,'YLabel'),'Rotation',0);
  else
    set(gca,'YTick', tiYs,'Yticklabel',[],'YDir','normal', 'FontSize',10)
  end

  if panel > (row-1)*col
    set(gca,'XTick', tiXs, 'YDir','normal', 'FontSize',10)
    xlabel('P_1(0)')
  else
    set(gca,'XTick', tiXs,'Xticklabel',[], 'YDir','normal')
  end

  set(gca,'YLim', [0, 1]);
  caxis([0, 3])
  
  % Panel titles
%  if isfinite(g)
%    title(strcat('\nu = ', num2str(g), '; c = ', num2str(c(1))), 'FontSize',10)
%  else
%    title(strcat('\nu = \infty', '; c = ', num2str(c(1))), 'FontSize',10)
%  end
  if isfinite(g)
    title(strcat('(', alphabet(panel), ': \nu=', num2str(g), ', c=', num2str(c(1)), ')'), 'FontSize',10)
  else
    title(strcat('(', alphabet(panel), ': \nu=\infty', ', c=', num2str(c(1)), ')'), 'FontSize',10)
  end  

  % Statistics
  final_mean=[mean2(Pop1); mean2(Pop2); mean2(Pop3); mean2(Pref)];
  final_stdv=[std2(Pop1);  std2(Pop2);  std2(Pop3);  std2(Pref)];
  %text(10,30,strcat(equil_str, num2str(final_mean,'%.2f'),'\pm', num2str(final_stdv,'%.2f')),'FontSize',5)

end

print -depsc coex_init_plants_pref.eps
print -dpng coex_init_plants_pref.png
