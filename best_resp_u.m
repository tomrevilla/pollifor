function u = best_resp_u(x, lambda, b, w)

% This function calculates the best response solution 
% for the consumpẗion preferences in flexible pollinators
% u: fractional preference for plant 1
% x: vector of the two plant and the animal densities
% lambda: vector of maximum benefits
% b: vector of consumption (e.g. pollination, frugivory) rates
% w: vector of resource loss rates

u = (lambda(1)*b(1)*x(1)*(w(2)+b(2)*x(3))-lambda(2)*b(2)*x(2)*w(1))/ ...
                  (b(1)*b(2)*x(3)*(lambda(1)*x(1) + lambda(2)*x(2)));

% Correction to keep fraction between 0 and 1
u = max(0, min(1, u));

% Correction to initialize replication equation
u = max(0.001, min(0.999, u));

end