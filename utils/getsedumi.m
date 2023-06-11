function [A, b, c, K] = getsedumi(X, ptype)
% Converts the optimal diagonal pre-conditioning problem into Sedumi format

if nargin == 1
    ptype = "L";
end % End if

X = X / trace(X);

if ptype == "L"
    
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
    
elseif ptype == "M"
    
    [n, ~] = size(X);
    
    K.s = [n, n];
    b = zeros(3, 1); b(1) = 1;
    
    A0 = [sparse(n * n, 1); vec(X)]';
    
    Ddiag = diag(X);
    [Druiz, Dtmp] = ruizscale(X, 100); %#ok
    Druiz = diag(Druiz);    
    
%     getoptcombprecond(X, diag(Ddiag), diag(Druiz));
    
    i = [ones(1, 2 * n), ones(1, 2 * n) * 2];
    idx = 1:n;
    j = [idx * n + idx - n, idx * n + idx - n + n * n, idx * n + idx - n, idx * n + idx - n + n * n];
    v = [Ddiag; -Ddiag; Druiz; -Druiz];
    
    A = sparse(i, j, v, 2, 2 * n * n, 4 * n);
    A = [A0; A];
    
    c = zeros(2 * n * n, 1);
    c(1:n * n) = vec(X);
    
    writesdpa('test.dat-s', A, b, c, K);
    
elseif ptype == "R"
    
    X = full(X);
    [m, n] = size(X);
    K.s = [n, n];
    b = zeros(m + 1, 1); b(1) = 1;
    A = zeros(m + 1, 2 * n * n);
    veye = vec(eye(n));
    A(1, 1:n * n) = veye;
    
    for i = 1:m
        xi = X(i, :);
        xixiT = vec(xi' * xi)';
        A(i + 1, 1:n * n) = -xixiT;
        A(i + 1, n * n + 1:end) = xixiT;
    end % End for
    
    A = sparse(A);
    c = zeros(2 * n * n, 1);
    c(n * n + 1 : end) = veye;
    
else
    error("Unsupported pre-conditioner type.");
end % End if

end % End function