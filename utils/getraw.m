function [data] = getraw(X, perturb)
% Equilibrate and get raw data from X 

XTX = X' * X;
data.X = X;
data.M = XTX;

data.perturb = 0.0;
data.cond = condest(XTX);

% Make the matrix well-defined and possible to do pre-conditioning
% We keep adding diagonal perturbations until the matrix is positive
% definite and the condition number <= 1e+08
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