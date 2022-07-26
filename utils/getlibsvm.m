function [] = getlibsvm(fname, perturb)

target = fullfile("datasets", "libsvm");

if isfile(fullfile(target, "mat", fname + "-opt.mat"))
    d = load(fullfile(target, "mat", fname + "-opt.mat"));
    data = d.data;
    fprintf("%30s %8d %10.3e %10.3e \n", fname, size(data.M, 1), data.cond, data.perturb);
    return;
end % End if

[~, X] = libsvmread(char(fullfile(target, 'meta', fname)));
data = getraw(X, perturb);
data.name = fname;

save(fullfile(target, "mat", fname + "-opt.mat"), "data");

% Build SDPA LHS
[A, b, c, K] = getsedumi(data.M, "L");
param.printlevel = 0;
writesdpa(fullfile(target, "sdp", fname + "-L.dat-s"), A, b, c, K, param);

% Build SDPA RHS
[A, b, c, K] = getsedumi(data.X, "R");
writesdpa(fullfile(target, "sdp", fname + "-R.dat-s"), A, b, c, K, param);

% matrix  dim  cond  perturb
fprintf("%30s %8d %10.3e %10.3e \n", fname, size(data.M, 1), data.cond, data.perturb);

end % End function