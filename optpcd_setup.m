% Set up the matlab toolbox for optimal preconditioner computation
% 
%       min_{e} cond(XE' * XE)
%       min_{d} cond(DX' * DX)
%       min_{d, e} cond(DXE' * DXE)
% 
% Paper available at https://arxiv.org/abs/2209.00809
% 
% Check README.md for more details

fprintf("Setting up the optimal diagonal preconditioning repo \n");
fprintf("The repo depends on:  \n\n");
fprintf("1. CVX   matlab toolbox    at http://cvxr.com/cvx/\n");
fprintf("2. HDSDP binary (optional) at https://github.com/COPT-Public/HDSDP \n\n");
fprintf("To reproduce the experiments, execute \n\n");
fprintf("    git checkout opt-precond \n\n");
fprintf("and ensure that the following data repos are installed \n\n");
fprintf("3. LIBSVM datasets              at https://www.csie.ntu.edu.tw/~cjlin/libsvm/\n");
fprintf("4. SuiteSparse matlab interface at https://sparse.tamu.edu/interfaces\n");
addpath(genpath("."));

fprintf("\nRunning a toy example.\n\n");
try
    opt_example
catch
    fprintf("Installation fails. :( \n");
    fprintf("Please check if cvx is available and installed correctly\n");
end % End try

fprintf("\nInstallation completes. Check README.md for usage details. \n");