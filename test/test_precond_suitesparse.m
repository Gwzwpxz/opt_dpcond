% Experiment on optimal preconditioning for suite-sparse datasets
clear; clc; close all;

ss_all = ssget;
max_col = 500;

% Get matrices with ncols <= max_col
ids = find(ss_all.ncols <= max_col);

fprintf("%30s %8s %10s %10s \n", "Matrix", "Size", "Cond.", "Perturb");
for i = 1:length(ids)
    getsuitesparse(ids(i), ss_all.Name{ids(i)});
end % End for

bad_list = [137];
fprintf("%30s %6s %6s %6s %s \n", "Mat", "Size", "cbef", "caft", "reduce");
for i = 1:length(ids)
    name = ss_all.Name{ids(i)};
    if ismember(i, bad_list)
        continue;
    end % End if
    d = load(fullfile("datasets", "suitesparse", "mat", name + "-opt.mat"));
    testoptprecond(d.data, "suitesparse", "hdsdp");
end % End for

