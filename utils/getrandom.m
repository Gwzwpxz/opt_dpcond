function [] = getrandom(n, sps, perturb)
% This function generates a test instance for optimal diagonal
% pre-conditioning from a randomly generated matrix

if nargin == 1
    perturb = 1e-06;
end % End if

target = fullfile("datasets", "random");
matname = sprintf("diag-bench-%d-%3.3e", n, sps);


if isfile(fullfile(target, "mat", matname + "-opt.mat"))
    d = load(fullfile(target, "mat", matname + "-opt.mat"));
    data = d.data;
    fprintf("%60s %8d %10.3e %10.3e \n", matname, size(data.M, 1), data.cond, data.perturb);
    return;
end % End if

rng(24);
X = sprandn(n, n, sps);

data = getraw(X, perturb);
data.name = matname;
save(fullfile(target, "mat", matname + "-opt.mat"), "data");

% Build SDPA LHS
[A, b, c, K] = getsedumi(data.M, "L");
param.printlevel = 0;
writesdpa(fullfile(target, "sdp", matname + "-L.dat-s"), A, b, c, K, param);

% Build SDPA RHS
[A, b, c, K] = getsedumi(data.X, "R");
writesdpa(fullfile(target, "sdp", matname + "-R.dat-s"), A, b, c, K, param);

% matrix  dim  cond  perturb
fprintf("%60s %8d %10.3e %10.3e \n", matname, size(data.M, 1), data.cond, data.perturb);

end % End function