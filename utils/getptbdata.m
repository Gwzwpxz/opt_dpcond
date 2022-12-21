function [data] = getptbdata(X, perturb)
% Perturb X' * X till its condition number is reasonable

XTX = X' * X;
data.X = X;
data.M = XTX;

data.perturb = 0.0;
data.cond = condest(XTX);

% Make the matrix better-conditioned
ptb = 0.0; count = 0;
while data.cond > 1e+06
    ptb = ptb + perturb;
    data.M = XTX + ptb * speye(size(XTX, 1));
    data.cond = condest(data.M);
    count = count + 1;
    if count > 3
        perturb = perturb * 2;
    end % End if
end % End if

data.perturb = ptb;
end % End function