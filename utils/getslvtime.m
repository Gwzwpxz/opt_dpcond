function [tpcg, tdirect, itcg, itpcg] = getslvtime(A, PA, ntest)

n = size(A, 1);
maxiter = 10 * n;

warning off;

itcg = 0;
itpcg = 0;

rng(100);

% Time sparse direct method
tic;
[L, ~, pmt] = lchol(A);
LT = L';
pinv = 1:n; pinv(pmt) = 1:n;
rhs = randn(n, 1);
for i = 1:ntest
    b = rhs(pmt, :);
    x = L \ b;
    x = LT \ x;
    x = x(pinv);
end % End for
tdirect = toc;

% A = full(A);
%
% tic;
% L = chol(A)';
% LT = L';
% rng(24);
% rhs = randn(n, 1);
% for i = 1:ntest
%     x = L \ rhs;
%     x = LT \ x;
% end % End for
% tdirect = toc;

tol = 1e-04; % max(norm(A * x - rhs), 1e-04);

% if length(nonzeros(PA)) > 0.7 * n * n
%     PA = full(PA);
% end % End if
% Time iterative method
PA = sparse(PA);
rhs10 = repmat(rhs, 1, 128);
tic;
for i = 1:(ntest / 128)
    [x] = rci(PA, rhs10, tol / 100, maxiter);
end % End for
tpcg = toc;

if norm(PA * x(1:n, 1) - rhs) > tol * norm(rhs)
    tpcg = 9.9;
else
    %     keyboard;
end % End if

end % End function