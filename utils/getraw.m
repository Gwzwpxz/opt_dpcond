function [data] = getraw(X, perturb)
% Equilibrate and get raw data from X 

XTX = X' * X;
data.X = X; % + sprandsym(11107, 0.00000005) * 1e+04;
data.M = XTX;

data.perturb = 0.0;
data.cond = condest(XTX);

% Make the matrix more well-defined
ptb = 0.0; count = 0;
while data.cond > 1e+08
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