function [D1, D2, feas] = d1d2SDP(X, kappa)

[m, n] = size(X);

cvx_begin sdp quiet
cvx_solver sdpt3
variable D1(n, n) diagonal
variable D2(m, m) diagonal
minimize 0
subject to
X' * D2 * X >= D1;
kappa * D1 >= X' * D2 * X;
D1 - speye(n) >= 0;
cvx_end % End cvx

feas = true;

if cvx_status ~= "Solved"
    feas = false;
end % End if

end % End function