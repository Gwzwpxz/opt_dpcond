% Experiment on optimal preconditioning for suite-sparse datasets
clear; clc; close all;

ss_all = ssget;
max_col = 1000;
max_row = 1000;
max_nnz = 25000;

% Get matrices with ncols <= max_col
idxr = find(ss_all.nrows <= max_row);
ids = find(ss_all.ncols >= 0);
idxc = find(ss_all.ncols <= max_col);
idxsp = find(ss_all.nnz < max_nnz);
ids = intersect(ids, idxr);
ids = intersect(ids, idxc);
ids = intersect(ids, idxsp);

fprintf("%30s %8s %10s %10s \n", "Matrix", "Size", "Cond.", "Perturb");
for i = 1:length(ids)
    getsuitesparse(ids(i), ss_all.Name{ids(i)});
end % End for

testsolve = true;
for i = 1:length(ids)
    name = ss_all.Name{ids(i)};
    d = load(fullfile("datasets", "suitesparse", "mat", name + "-opt.mat"));
    testoptprecond(d.data, "suitesparse", "cvx", "M", testsolve);
end % End for


% bad_list = ["qc324"];
% fileID = fopen("sslog-LR.txt", 'w+');
% fprintf(fileID, "%4s %30s %6s %6s %6s %6s %6s %s\n", "Type", "Mat", "Size", "cbef", "caft", "cruiz", "reduce", "time");
% for i = 1:length(ids)
%     name = ss_all.Name{ids(i)};
%     fprintf("Running %s \n", name);
%     if ismember(name, "qc324")
%         continue;
%     end % End if
%     d = load(fullfile("datasets", "suitesparse", "mat", name + "-opt.mat"));
%     testoptprecond(d.data, "suitesparse", "hdsdp", "LR", fileID);
% end % End for
% fclose(fileID);

% bad_list = [];
% fileID = fopen("sslog-L-normalized.txt", 'w+');
% fprintf(fileID, "%4s %30s %6s %6s %6s %6s %s\n", "Type", "Mat", "Size", "cbef", "caft", "reduce", "time");
% for i = 1:length(ids)
%     name = ss_all.Name{ids(i)};
%     if ismember(name, "qc324")
%         continue;
%     end % End if
%     d = load(fullfile("datasets", "suitesparse", "mat", name + "-opt.mat"));
%     testoptprecond(d.data, "suitesparse", "cvx", "L");
% end % End for
% fclose(fileID);

% ids = intersect(ids, idxr);
% bad_list = [224];
% fileID = fopen("sslog-R.txt", 'w+');
% fprintf(fileID, "%4s %30s %6s %6s %6s %6s %s\n", "Type", "Mat", "Size", "cbef", "caft", "reduce", "time");
% for i = 1:length(ids)
%     name = ss_all.Name{ids(i)};
%     if ismember(name, "qc324")
%         continue;
%     end % End if
%     d = load(fullfile("datasets", "suitesparse", "mat", name + "-opt.mat"));
%     testoptprecond(d.data, "suitesparse", "hdsdp", "R", fileID);
% end % End for
% fclose(fileID);
