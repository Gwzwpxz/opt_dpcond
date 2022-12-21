function [D1, D2] = gettwosidedprecond(X, ub)
% Compute two-sided pre-conditioner using SDP
if nargin == 1
    ub = 1e+06;
end % End if

lb = 1;
feas = false;

while true
    [D1, D2, feas] = iskappafeas(X, ub);
    if feas || ub > 1e+10
        break;
    end % End if
    ub = ub * 10;
end % End while

fprintf("Solving a two-sided preconditioning problem using bisection \n");
fprintf("%10s %10s \n", "kappa", "ub - lb");
while ub - lb > 1e-02 * (abs(ub) + abs(lb) + 1)
    diff = ub - lb;
    kappa = lb + diff / 2;
    [D1, D2, feas] = iskappafeas(X, kappa);
    if feas
        ub = kappa;
    else
        lb = kappa;
    end
    fprintf("%10.3e %10.3e \n", kappa, diff);
end % End while

if ~feas
    [D1, D2, ~] = iskappafeas(X, ub);
end % End if

end % End function