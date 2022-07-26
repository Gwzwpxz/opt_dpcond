% Experiment on optimal preconditioning for suite-sparse datasets
clear; clc; close all;

ss_all = ssget;
max_col = 100;
max_row = 100;

% Get matrices with ncols <= max_col
idxr = find(ss_all.nrows > max_row);
ids = find(ss_all.ncols > max_col);
idxc = find(ss_all.ncols <= 1000);
ids = intersect(ids, idxr);
ids = intersect(ids, idxc);

fprintf("%30s %8s %10s %10s \n", "Matrix", "Size", "Cond.", "Perturb");
for i = 1:length(ids)
    getsuitesparse(ids(i), ss_all.Name{ids(i)});
end % End for

bad_list = ["qc324"];
fileID = fopen("sslog-LR.txt", 'w+');
fprintf(fileID, "%4s %30s %6s %6s %6s %6s %6s %s\n", "Type", "Mat", "Size", "cbef", "caft", "cruiz", "reduce", "time");
for i = 1:length(ids)
    name = ss_all.Name{ids(i)};
    fprintf("Running %s \n", name);
    if ismember(name, "qc324")
        continue;
    end % End if
    d = load(fullfile("datasets", "suitesparse", "mat", name + "-opt.mat"));
    testoptprecond(d.data, "suitesparse", "hdsdp", "LR", fileID);
end % End for
fclose(fileID);

% bad_list = ["qc324"];
% fileID = fopen("sslog-L.txt", 'w+');
% fprintf(fileID, "%4s %30s %6s %6s %6s %6s %s\n", "Type", "Mat", "Size", "cbef", "caft", "reduce", "time");
% for i = 1:length(ids)
%     name = ss_all.Name{ids(i)};
%     if ismember(name, "qc324")
%         continue;
%     end % End if
%     d = load(fullfile("datasets", "suitesparse", "mat", name + "-opt.mat"));
%     testoptprecond(d.data, "suitesparse", "hdsdp", "L", fileID);
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
