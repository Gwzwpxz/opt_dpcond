clear;
clc;
close all;

% rng(1);

% a = 1e+05;
% c = 1;
% b = sqrt(a * c - 1);
% A = [a, b; b, c];
% A = A / trace(A);

% v = (trace(A) / 2 + sqrt(trace(A)^2 / 4 - (a * c - b^2)))^2 / det(A);
% 
% v2 = ((sqrt(a * c) + b)) / (sqrt(a * c) - b);
% 
% A2 = [1, b / sqrt(a * c);
%       b / sqrt(a * c), 1];
% 
% [Dopt] = getcvxdiag(A, 'L');
% 
% Aopt = sqrt(Dopt) \ (sqrt(Dopt) \ full(A))';

A = genmat(200);

Djacob = diag(diag(A));
% [Druiz, ~] = ruizscale(A, 100);
% Druiz = diag(1./diag(Druiz.^2));
% Dspinv = diag(diag(A) ./ sum(A.^2, 2));
% Dspinv = diag(1./diag(Dspinv));

Ajacob = sqrt(Djacob) \ (sqrt(Djacob) \ full(A))';
% Aruiz = sqrt(Druiz) \ (sqrt(Druiz) \ full(A))';
% Aspinv = sqrt(Dspinv) \ (sqrt(Dspinv) \ full(A))';

% cond_spinv = cond(full(Aspinv));
% cond_ruiz = cond(full(Aruiz));

[D, cvx_time] = getoptcombprecond(A, Djacob, eye(size(A, 1)));
[Dopt] = getcvxdiag(A, 'L');

Aopt = sqrt(Dopt) \ (sqrt(Dopt) \ full(A))';
Acomb = sqrt(D) \ (sqrt(D) \ full(A))';

cond_jacob = cond(full(Ajacob));
cond_opt = cond(full(Aopt));
cond_comb = cond(full(Acomb));

fprintf("Jacobi: %e Opt: %e Combined: %e \n", cond_jacob, cond_opt, cond_comb);

if cond_jacob > 100 * cond_opt
    keyboard;
end % End if 

rhs = randn(size(A, 1), 1);
pcg(Aopt, rhs, 1e-10, 10000);
pcg(Acomb, rhs, 1e-10, 10000);
pcg(Ajacob, rhs, 1e-10, 10000);

function [A] = genmat(n)

i = 1:n;
j = 1:n;
v = max(logspace(-4, 1, n)', 1e-08);
A = sparse(i, j, v(randperm(n)));

a = randn(n, 1) * 20;

A = A + a * a';

% A = A + sprandsym(n, n, 0.1) * 100;
% A(1:n, 1:n) = max(A(1:n, 1:n), 1e-04);

% A = A' * A;
% while 1
%     try 
%         chol(A);
%         break;
%     catch
%         A = A + eye(n);
%     end % End try
% end % End while

A = A / trace(A);

end % End function

