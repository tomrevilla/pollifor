% REPLICATOR MODEL SIMULATION SCRIPT -- SCENARIO 2
%
% x: vector of initial population densities and preference
% r: vector of plant yields
% e: vector of animal yields
% b: vector of consumption (pollination or frugivory) rates
% a: vector of (flower, fruit) supply rates
% w: vector of loss rates
% m: vector of plant and animal mortality rates
% c: vector of plant competition coefficients
% K: vector of plant zero yield density
% g: foraging adaptation rate
% pref: preference model cases (see odesystem.m)

tic

% Default parameters
r = 0.1*[1, 1];
e = 0.1*[2, 1];
b = 0.1*[1, 1];
a = 0.4*[1, 1];
w = 0.25*[1, 1];
m = [0.01, 0.01, 0.1];
c = 0.4*[1, 1];         % competition level
K = 50*[1, 1];
g = 0.1;                % adaptation level

% Operation parameters
% pref = 0: use with g = 0,   fixed preference
%        1: use with g = inf, ESS preference
%        3: use with g > 0,   replicator equation
%        4: use with g > 0,   best response dynamics (not employed)

pref = 3;               % calls 'odesystem' with replicator eq. option
gridres = 100;          % size of the initial condition grid
extinc_level = 1e-6;
options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6 1e-6], 'MaxStep', 1e-10);
tspan = 20000;

% Plant population density grid: The sum of plant densities is constrained to be 
% equal to the average of their carrying capacities, so we run across plant 1 
% density. The preference runs from 0 to 1. Pollinator density is A=2
Pop1 = zeros(gridres);
Pop2 = zeros(gridres);
Pop3 = zeros(gridres);
Pref = zeros(gridres);
init_pla = linspace(0, ceil(mean(K)), gridres);
init_pol = 2;
init_prf = linspace(0, 1, gridres);
[X,Y] = meshgrid(init_pla, init_prf);

% Parallel runs. The variable 'pref' determines the preference model to run.
matlabpool open local 7

parfor xx = 1:gridres
  dummy = zeros(4,gridres);
  for yy = 1:gridres
        n0 = [init_pla(xx), K(1)-init_pla(xx), init_pol, init_prf(yy)];
        [t, y] = ode45(@odesystem, [0, tspan], n0, [], ...
                 r, e, b, a, w, m, c, K, g, pref);
    
    disp([n0; y(end,:)])   % see progress ...

    dummy(1,yy) = y(end,1);
    dummy(2,yy) = y(end,2);
    dummy(3,yy) = y(end,3);
    dummy(4,yy) = y(end,4);
  end
  Pop1(xx,:) = dummy(1,:);
  Pop2(xx,:) = dummy(2,:);
  Pop3(xx,:) = dummy(3,:);
  Pref(xx,:) = dummy(4,:);
end

matlabpool close

elapsed_time = toc

% save all. Choose appropriate /path/name.mat. E.g., see 'datastring' definition 
% in 'coex_init_plants_pref.m'
save('data_scen2.mat')
clear
