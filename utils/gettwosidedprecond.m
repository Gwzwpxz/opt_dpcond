function [D1, D2] = gettwosidedprecond(X, ub)
% Compute two-sided pre-conditioner using SDP
if nargin == 1
    ub = 1e+06;
end % End if

% [refL, refR] = getitertwosidedprecond(X);

lb = 1;
feas = false;

while true
   [D1, D2, feas] = d1d2SDP(X, ub);
   if feas || ub > 1e+10
       break;
   end % End if
   ub = ub * 10;
end % End while

if ub >= 1e+10
  [m, n] = size(X);
  D1 = speye(n);
  D2 = speye(m);
end 

% fprintf("%10s %10s \n", "kappa", "ub - lb");
while ub - lb > 1e-05 * (abs(ub) + abs(lb) + 1)
    diff = ub - lb;
    kappa = lb + diff / 2;
    [D1, D2, feas] = d1d2SDP(X, kappa);
    if feas
       ub = kappa;
    else
       lb = kappa;
    end
%     fprintf("%10.3e %10.3e \n", kappa, diff);
end % End while

if ~feas 
    [D1, D2, ~] = d1d2SDP(X, ub);
end % End if

end % End function

