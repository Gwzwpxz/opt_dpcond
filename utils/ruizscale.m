function [D, E] = ruizscale(A, maxiter)
% Implement the Ruiz scaling algorithm for general matrix

if nargin == 1
    maxiter = 100;
end % End if

[m, n] = size(A);
D = ones(m, 1);
E = ones(n, 1);

for i = 1:maxiter
    
    dR = sqrt(sum(abs(A), 2));
    dC = sum(abs(A), 1).^(-1/2);
%     dR = sqrt(max(abs(A), [], 2));
%     dC = max(abs(A)).^(-1/2);
    R = diag(dR);
    C = diag(dC);
    A = R \  (A * C);
    D = D ./ dR;
    E = E .* dC';
    
    if norm(dR - 1, 'inf') <= 1e-10 && norm(dC - 1, 'inf') <= 1e-10
        break;
    end % End if
    
end % End for

D = sparse(1:m, 1:m, D);
E = sparse(1:n, 1:n, E);

end % End function