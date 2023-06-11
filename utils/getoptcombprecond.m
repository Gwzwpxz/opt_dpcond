function [D, cvx_time] = getoptcombprecond(X, D1, D2, D3)
% Get a pre-conditioner that is the linear/affine combination of two base
% pre-conditioners

D1 = D1 / trace(D1);
D2 = D2 / trace(D2);
X = X / trace(X);

if nargin == 4
    
    D3 = D3 / trace(D3);
    
    cvx_begin sdp quiet
    cvx_solver mosek
    variable beta1 % nonnegative
    variable beta2 % nonnegative
    variable beta3 % nonnegative
    variable tau
    maximize tau
    subject to
    % tau <= 1.0; %#ok
    X - beta1 * D1 - beta2 * D2 - beta3 * D3 >= 0; %#ok
    beta1 * D1 + beta2 * D2 + beta3 * D3 - X * tau >= 0; %#ok
    tic;
    cvx_end
    
    D = beta1 * D1 + beta2 * D2 + beta3 * D3;
    
else
    
    cvx_begin sdp quiet
    cvx_solver mosek
    variable beta1 % nonnegative
    variable beta2 % nonnegative
    variable tau
    maximize tau
    subject to
    X - beta1 * D1 - beta2 * D2 >= 0; %#ok
    beta1 * D1 + beta2 * D2 - X * tau >= 0; %#ok
    tic;
    cvx_end
    
    D = beta1 * D1 + beta2 * D2;
    
end % End if
cvx_time = toc;

% fprintf("%f seconds\n", cvx_time);

end % End function