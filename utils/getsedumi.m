function [A, b, c, K] = getsedumi(X)
% Converts the optimal diagonal pre-conditioning problem into Sedumi format

[n, ~] = size(X);

% Two SDP blocks
K.s = [n, n]; 
% Dual obj, b(1) is kappa 
b = zeros(n + 1, 1); b(1) = 1;
% Conic constriants

A0 = [sparse(n * n, 1); vec(X)]';

i = [1:n, 1:n];
idx = 1:n;
j = [idx * n + idx - n, idx * n + idx - n + n * n]; 
v = [ones(n, 1); ones(n, 1) * -1];
A = sparse(i, j, v, n, 2 * n * n, 2 * n);

A = [A0; A];

c = zeros(2 * n * n, 1);
c(1:n * n) = vec(X);

end % End function