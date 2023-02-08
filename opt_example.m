clear; close all;

% Load data
load(fullfile('.', 'data', 'diag-bench-100-1.000e-01-opt.mat'));
X = data.X;

% Setup problem
param.ptype = 'R';
param.perturb = true; 
Eprob = getoptprob(X, param);

param.ptype = 'L';
Dprob = getoptprob(X, param);

param.ptype = 'T';
DEprob = getoptprob(X, param);

% Solve problem
Eopt = optprecond(Eprob);
Dopt = optprecond(Dprob);
DEopt = optprecond(DEprob);

% Get performance
fprintf("Condition number of X  : %5.2e \n", cond(full(X)));
fprintf("Condition number of XE : %5.2e \n", cond(full(Eopt.pX)));
fprintf("Condition number of DX : %5.2e \n", cond(full(Dopt.pX)));
fprintf("Condition number of DXE: %5.2e \n", cond(full(DEopt.pX)));

% Export dat-s
exportoptprob(Eprob, fullfile('.'), 'Eprob');