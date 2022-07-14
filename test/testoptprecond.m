function [D] = testoptprecond(data, type, solver, fileID)

if nargin == 3
    fileID = [];
end % End if

DSDP_RETCODE_OK = 2;

M = data.M;
cond_before = cond(full(M), 2);
sdptime = 0.0;

if solver == "cvx"
    D = getcvxdiag(M);
else
    target = fullfile("datasets", type, "sdp");
    hdsdppath = fullfile("hdsdp", "optprecond");
    targetsdp = fullfile(target, data.name + ".dat-s");
    targetcsv = fullfile(target, data.name + ".csv");
    if ~isfile(targetcsv)
        tic;
        cmd = sprintf("%s %s", hdsdppath, targetsdp);
        retcode = system(cmd);
        if retcode ~= DSDP_RETCODE_OK
            return;
        end % End if
        sdptime = toc;
    end % End if
    y = csvread(targetcsv);
    D = diag(y(2:end - 1));
end % End if

cond_after = cond(sqrt(D) \ (sqrt(D) \ full(data.M))');
reduce_cond = cond_before - cond_after;

if reduce_cond < 0
    return;
end % End if

rdc = reduce_cond / cond_before;
class = 3;
if rdc >= 0.8
    class = 1;
elseif rdc >= 0.5
    class = 2;
end % End if

log = sprintf("%d %30s %6d %6.3e %6.3e %f %f\n", class, data.name,...
    size(D, 1), cond_before, cond_after, reduce_cond / cond_before, sdptime);

if isempty(fileID)
    fprintf(log);
else
    fprintf(fileID, log);
end % End if

end % End function
