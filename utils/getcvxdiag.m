function [D] = getcvxdiag(X)
% Get the optimal diagonal pre-conditioner using CVX

[n, ~] = size(X);
cvx_begin sdp
cvx_solver sdpt3
variable d(n, 1) nonnegative
variable tau nonnegative
maximize tau
subject to
    tau <= 1.0;
    X - diag(d) >= 0;
    diag(d) - X * tau >= 0;
cvx_end

% cvx_begin sdp
% cvx_solver sdpt3
% variable d(n, 1) nonnegative
% variable tau
% maximize tau
% subject to
%     X - diag(d) >= 0;
%     diag(d) - X * tau >= 0;
% cvx_end

D = diag(d);

end % End function