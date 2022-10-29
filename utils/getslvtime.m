function [tpcg, tdirect, itcg, itpcg] = getslvtime(A, PA, ntest)

n = size(A, 1);
maxiter = 10000;

warning off;

% Time sparse direct method
% tic; 
% [L, ~, pmt] = lchol(A);
% LT = L';
% pinv = 1:n; pinv(pmt) = 1:n;
% rhs = randn(n, 1);
% for i = 1:ntest
%     b = rhs(pmt, :);
%     x = L \ b;
%     x = LT \ x;
%     x = x(pinv);
% end % End for
% tdirect = toc;

A = full(A);

tic; 
L = chol(A)';
LT = L';
rng(24);
rhs = randn(n, 1);
for i = 1:ntest
    x = L \ rhs;
    x = LT \ x;
end % End for
tdirect = toc;

tol = 1e-05; % max(norm(A * x - rhs), 1e-04);

% if length(nonzeros(PA)) > 0.7 * n * n
%     PA = full(PA);
% end % End if
% Time iterative method
tic;
for i = 1:ntest
[x, flag, r, itpcg] = pcg(PA, rhs, tol, maxiter);
end % End for
tpcg = toc;

[x, flag, ~, itcg] = pcg(A, rhs, tol, maxiter);

if itpcg > itcg
%     keyboard;
end % End if

end % End function