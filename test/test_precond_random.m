% Experiment on optimal preconditioning for suite-sparse datasets
clear; clc; close all;

nsize = [10, 50, 100, 500];
sps = [0.8, 0.6, 0.5, 0.3];

% Generate random instances
fprintf("%30s %8s %10s %10s \n", "Matrix", "Size", "Cond.", "Perturb");
for i = 1:length(nsize)
    n = nsize(i);
    getrandom(n, sps(i), 1e-06);
end % End for

% fprintf("%30s %6s %6s %6s %s \n", "Mat", "Size", "cbef", "caft", "reduce");
% for i = 1:length(nsize)
%     name = sprintf("diag-bench-%d-%3.3e", nsize(i), sps(i));
%     d = load(fullfile("datasets", "random", "mat", name + "-opt.mat"));
%     testoptprecond(d.data, "random", "hdsdp", "L");
% end % End for
% 
% fprintf("%30s %6s %6s %6s %s \n", "Mat", "Size", "cbef", "caft", "reduce");
% for i = 1:length(nsize)
%     name = sprintf("diag-bench-%d-%3.3e", nsize(i), sps(i));
%     d = load(fullfile("datasets", "random", "mat", name + "-opt.mat"));
%     testoptprecond(d.data, "random", "hdsdp", "R");
% end % End for
% 
for i = 1:length(nsize)
    name = sprintf("diag-bench-%d-%3.3e", nsize(i), sps(i));
    d = load(fullfile("datasets", "random", "mat", name + "-opt.mat"));
    testoptprecond(d.data, "random", "cvx", "M");
end % End for






