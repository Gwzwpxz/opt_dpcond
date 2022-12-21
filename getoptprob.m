function [prob] = getoptprob(X, params)
% Build a problem structure describing a diagonal preconditioning problem
% X      : m by n matrix
% Parameters:
% ptype: 'L', 'R', 'T' for left (DX), right (XE) and two-sided (DXE)
% perturb: swictch of perturbation of diagonal of X'X when right 
%          preconditioning E(X'X)E is used. Otherwise ignored.
%          If taken < 0, then perturbation will continue till 
%          cond(X'X + eps * I) < 1e+06
if nargin == 1
    params.ptype = 'R';
    params.perturb = true;
end % End if

ptype = upper(params.ptype);

if ptype == 'R'
    if params.perturb
        prob = getptbdata(X, 1.0);
    else
        prob.ptype = ptype;
        prob.cond = condest(X' * X);
        prob.X = X;
        prob.E = [];
    end % End if
    prob.E = [];
elseif ptype == 'L'
    prob.ptype = ptype;
    prob.cond = condest(X' * X);
    prob.X = X;
    prob.D = [];
else
    prob.ptype = 'T';
    prob.cond = condest(X' * X);
    prob.X = X;
    prob.D = [];
    prob.E = [];
end % End if

prob.ptype = ptype;
prob.optcond = inf;

end % End function