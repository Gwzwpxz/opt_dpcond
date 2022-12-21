function [D1, D2, feas] = iskappafeas(X, kappa) %#ok
% Check feasibility of an SDP bisection subproblem

[m, n] = size(X); %#ok

cvx_begin sdp quiet
% cvx_solver sdpt3
variable D1(n, n) diagonal
variable D2(m, m) diagonal semidefinite
minimize 0
subject to
    X' * D2 * X >= D1; %#ok
    kappa * D1 >= X' * D2 * X; %#ok
    D1 - speye(n) >= 0; %#ok
cvx_end % End cvx

feas = true;

if cvx_status ~= "Solved"
    feas = false;
end % End if

end % End function