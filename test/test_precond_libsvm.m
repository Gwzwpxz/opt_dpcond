% Experiment on optimal preconditioning for suite-sparse datasets
clear; clc; close all;

datasets = dir(fullfile("datasets", "libsvm", "meta"));
fnames = {datasets.name}';

fprintf("%30s %8s %10s %10s \n", "Matrix", "Size", "Cond.", "Perturb");
for i = 1:length(fnames)
    if ismember(fnames{i}, ['.', '..', '.DS_Store'])
        continue;
    end % End if
    getlibsvm(fnames{i}, 1e-06);
end % End for

fprintf("%30s %6s %6s %6s %s \n", "Mat", "Size", "cbef", "caft", "reduce");
for i = 1:length(fnames)
    if ismember(fnames{i}, ['.', '..', '.DS_Store'])
        continue;
    end % End if
    name = fnames{i};
    d = load(fullfile("datasets", "libsvm", "mat", name + "-opt.mat"));
    testoptprecond(d.data, "libsvm", "hdsdp");
end % End for
