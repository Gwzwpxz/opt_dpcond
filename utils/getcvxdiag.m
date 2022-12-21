function [D] = getcvxdiag(X, type)
% Get the optimal diagonal pre-conditioner using CVX

if nargin == 1
    type = "L";
end % End if

[m, n] = size(X);
if type == "L"
    cvx_begin sdp quiet
        cvx_solver sdpt3
        variable d(n, 1) nonnegative
        variable tau nonnegative
        maximize tau
        subject to
        tau <= 1.0;
        X - diag(d) >= 0;
        diag(d) - X * tau >= 0;
        d >= 0;
    cvx_end
elseif type == "R"
    cvx_begin sdp quiet
    cvx_solver sedumi
    variable d(m, 1) nonnegative
    variable tau
    maximize tau
    subject to
    tau * speye(n) - X' * diag(d) * X <= 0;
    X' * diag(d) * X  - speye(n) <= 0;
    d >= 0;
    cvx_end
else
    error("Invalid pre-conditioner type");
end % End if

D = diag(d);
% cvx_begin sdp
% cvx_solver sdpt3
% variable d(n, 1) nonnegative
% variable tau
% maximize tau
% subject to
%     X - diag(d) >= 0;
%     diag(d) - X * tau >= 0;
% cvx_end

end % End function