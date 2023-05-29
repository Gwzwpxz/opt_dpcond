function [D] = getcvxdiag(X, type)
% Get the optimal diagonal pre-conditioner using CVX

if nargin == 1
    type = "L";
end % End if

[~, n] = size(X);
if type == "L"
    cvx_begin sdp quiet
        cvx_solver sdpt3
        variable d(n, 1) nonnegative
        variable tau nonnegative
        maximize tau
        subject to
        tau <= 1.0; %#ok
        X - diag(d) >= 0; %#ok
        diag(d) - X * tau >= 0; %#ok
        d >= 0; %#ok
    cvx_end
elseif type == "R"
    cvx_begin sdp quiet
    cvx_solver sedumi
    variable d(m, 1) nonnegative
    variable tau
    maximize tau
    subject to
    tau * speye(n) - X' * diag(d) * X <= 0; %#ok
    X' * diag(d) * X  - speye(n) <= 0; %#ok
    d >= 0; %#ok
    cvx_end
else
    error("Invalid pre-conditioner type");
end % End if

D = diag(d);

end % End function