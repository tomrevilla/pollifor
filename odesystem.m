function xdot = odesystem(t, x, r, e, b, a, w, m, c, K, g, pref)

% Differential equations for the obligatory mutualism between two plants and 
% one pollinator or frugivore.
%
% x: vector of population densities and preference (*)
% r: vector of plant yields
% e: vector of animal yields
% b: vector of consumption (pollination or frugivory) rates
% a: vector of (flower, fruit) supply rates
% w: vector of loss rates
% m: vector of plant and animal mortality rates
% c: vector of plant competition coefficients
% K: vector of plant zero yield density
% g: foraging adaptation rate (*)
% pref: preference model (*)
%
% (*): Depending on the preference model we will have 3 or 4 differential 
% equations. In '0' and '1' preferences are fixed or flexible but both are 
% instantaneous, so we only use 3 equations. Cases '2' and '3' uses a
% differential equation for the preference, thus we use 4 equations and we 
% need a fourth intial condition for vector 'x'.
%
% The best response is called externally using the 'best_resp_u' function.

switch pref
  case 0
  % Fixed preferences  
  xdot = zeros(3,1);

  xdot(1) = (r(1)*(1-(x(1)+c(2)*x(2))/K(1))* ...
             b(1)*a(1)*x(3)/(w(1) + b(1)*x(3)) -m(1))*x(1);
              
  xdot(2) = (r(2)*(1-(x(2)+c(1)*x(1))/K(2))* ...
             b(2)*a(2)*x(3)/(w(2) + b(2)*x(3)) -m(2))*x(2);
          
  xdot(3) = (e(1)*b(1)*a(1)*x(1)/(w(1) + b(1)*x(3)) + ...
             e(2)*b(2)*a(2)*x(2)/(w(2) + b(2)*x(3)) -m(3))*x(3);

  case 1
  % Flexible preferences equal to the best response
  u1 = best_resp_u(x, e.*a, b, w);
  u2 = 1 - u1;
  
  xdot = zeros(3,1);

  xdot(1) = (r(1)*(1-(x(1)+c(2)*x(2))/K(1)) * ... 
             u1*b(1)*a(1)*x(3)/(w(1) + u1*b(1)*x(3)) -m(1))*x(1);
              
  xdot(2) = (r(2)*(1-(x(2)+c(1)*x(1))/K(2)) * ... 
             u2*b(2)*a(2)*x(3)/(w(2) + u2*b(2)*x(3)) -m(2))*x(2);
          
  xdot(3) = (e(1)*u1*b(1)*a(1)*x(1)/(w(1) + u1*b(1)*x(3)) + ... 
             e(2)*u2*b(2)*a(2)*x(2)/(w(2) + u2*b(2)*x(3)) -m(3))*x(3);
  
  case 2
  % Flexible preferences with best response dynamics (it was not used)
  xdot = zeros(4,1);

  xdot(1) = (r(1)*(1-(x(1)+c(2)*x(2))/K(1)) * ...
             x(4)*b(1)*a(1)*x(3)/(w(1) + x(4)*b(1)*x(3)) -m(1))*x(1);
              
  xdot(2) = (r(2)*(1-(x(2)+c(1)*x(1))/K(2)) *... 
             (1-x(4))*b(2)*a(2)*x(3)/(w(2) + (1-x(4))*b(2)*x(3)) -m(2))*x(2);
          
  xdot(3) = (e(1)*  x(4)  *b(1)*a(1)*x(1)/(w(1) +   x(4)  *b(1)*x(3)) + ...
             e(2)*(1-x(4))*b(2)*a(2)*x(2)/(w(2) + (1-x(4))*b(2)*x(3)) -m(3))*x(3);
  
  xdot(4) = g*(1/(1 + exp(100*(x(4) - best_resp_u(x, e.*a, b, w)))) - x(4));
  
  case 3
  % Flexible preferences using replicator equation
  xdot = zeros(4,1);

  xdot(1) = (r(1)*(1-(x(1)+c(2)*x(2))/K(1)) * ...
             x(4)*b(1)*a(1)*x(3)/(w(1) + x(4)*b(1)*x(3)) -m(1))*x(1);
              
  xdot(2) = (r(2)*(1-(x(2)+c(1)*x(1))/K(2)) *... 
             (1-x(4))*b(2)*a(2)*x(3)/(w(2) + (1-x(4))*b(2)*x(3)) -m(2))*x(2);
          
  xdot(3) = (e(1)*  x(4)  *b(1)*a(1)*x(1)/(w(1) +   x(4)  *b(1)*x(3)) + ...
             e(2)*(1-x(4))*b(2)*a(2)*x(2)/(w(2) + (1-x(4))*b(2)*x(3)) -m(3))*x(3);

  xdot(4) = g*x(4)*(1-x(4))*(e(1)*b(1)*a(1)*x(1)/(w(1) +   x(4)  *b(1)*x(3)) ...
                           - e(2)*b(2)*a(2)*x(2)/(w(2) + (1-x(4))*b(2)*x(3)));
                           
  
  otherwise
    disp('ERROR: no preference model was chosen.')

end


end
