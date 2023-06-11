function [] = getsuitesparse(id, matname, perturb)
% This function generates a test instance for optimal diagonal
% pre-conditioning from a suite-sparse matrix

if nargin == 2
    perturb = 1e-06;
end % End if

target = fullfile("datasets", "suitesparse");

if isfile(fullfile(target, "mat", matname + "-opt.mat"))
    p = load(fullfile(target, "mat", matname + "-opt.mat"));
    data = p.data;
    fprintf("%30s %8d %10.3e %10.3e \n", matname, size(data.M, 1), data.cond, data.perturb);
    return;
end % End if

% Process data and add perturbation
p = ssget(id);
X = p.A;
data = getraw(X, perturb);
data.name = matname;

save(fullfile(target, "mat", matname + "-opt.mat"), "data");

% Build SDPA for LHS pre-conditioning
[A, b, c, K] = getsedumi(data.M, "L");
param.printlevel = 0;
writesdpa(fullfile(target, "sdp", matname + "-L.dat-s"), A, b, c, K, param);

[A, b, c, K] = getsedumi(data.M, "M");
param.printlevel = 0;
writesdpa(fullfile(target, "sdp", matname + "-M.dat-s"), A, b, c, K, param);

% Build SDPA for RHS pre-conditioning
% [A, b, c, K] = getsedumi(data.X, "R");
% writesdpa(fullfile(target, "sdp", matname + "-R.dat-s"), A, b, c, K, param);

% matrix  dim  cond  perturb
fprintf("%30s %8d %10.3e %10.3e \n", matname, size(data.M, 1), data.cond, data.perturb);

end % End function