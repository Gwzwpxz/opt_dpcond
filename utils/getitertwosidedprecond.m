function [L, R] = getitertwosidedprecond(X)
% Iteratively compute the optimal preconditioning problem

for i = 1:20

[L] = getcvxdiag(X' * X, "L");
X = X * diag(diag(L).^(-0.5));

fprintf("%10.3e \n", cond(X' * X));
[R] = getcvxdiag(X, "R");
X = sqrt(R) * X;
fprintf("%10.3e \n", cond(X' * X));

end % End for



end % End function