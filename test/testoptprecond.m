function [D] = testoptprecond(data, type, solver, pcond, fileID)

if nargin == 4
    fileID = [];
end % End if

DSDP_RETCODE_OK = 2;
sdptime = 0.0;

if pcond == "L"
    M = data.M;
    cond_before = cond(full(M), 2);
    if solver == "cvx"
        D = getcvxdiag(M, pcond);
    else
        target = fullfile("datasets", type, "sdp");
        hdsdppath = fullfile("hdsdp", "optprecond");
        targetsdp = fullfile(target, data.name + "-L.dat-s");
        targetcsv = fullfile(target, data.name + "-L.csv");
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
    
elseif pcond == "R"
    % GET RHS pre-conditioning
    X = data.X; XTX = X' * X;
    cond_before = cond(full(XTX));
    if cond_before > 1e+08
        fprintf("Here, \n");
        return;
    else
        target = fullfile("datasets", type, "sdp");
        hdsdppath = fullfile("hdsdp", "optprecond");
        targetsdp = fullfile(target, data.name + "-R.dat-s");
        targetcsv = fullfile(target, data.name + "-R.csv");
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
        try
            D = diag(y(2:end - 1));
            cond_after = cond(full(X' * D * X));
        catch 
            return;
        end % End try
        reduce_cond = cond_before - cond_after;
    end % End if
    
elseif pcond == "LR"
    % Get two-sided pre-conditiong
    X = data.X; XTX = X' * X;
    cond_before = cond(full(XTX));
    if cond_before >= 1e+08
        return;
    end % End if
    [d1, d2] = gettwosidedprecond(X, 1e+04);
    D.d1 = d1;
    D.d2 = d2;
    XTX = X' * d2 * X;
    try
        cond_after = cond(sqrt(d1) \ (sqrt(d1) \ XTX)');
    catch
        return;
    end % End try
    reduce_cond = cond_before - cond_after;
end % End if

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

if pcond == "L"
    Ddiag = diag(diag(data.M));
    cond_diag = cond(sqrt(Ddiag) \ (sqrt(Ddiag) \ full(data.M))');
    log = sprintf("%d(%s) %30s %6d %6.3e %6.3e %6.3e %f %f\n", class, pcond, data.name,...
        size(D, 1), cond_before, cond_after, cond_diag, reduce_cond / cond_before, sdptime);
elseif pcond == "R"
    log = sprintf("%d(%s) %30s %6d %6.3e %6.3e %f %f\n", class, pcond, data.name,...
        size(D, 1), cond_before, cond_after, reduce_cond / cond_before, sdptime);
else
    [L, R] = ruizscale(data.X, 100);
    Xcond = L * data.X * R;
    cond_ruiz = cond(full(Xcond' * Xcond));
    log = sprintf("%d(%s) %30s %6d %6.3e %6.3e %6.3e %f\n", class, pcond, data.name,...
        size(X, 1), cond_before, cond_after, cond_ruiz, reduce_cond / cond_before);
end % End if

if isempty(fileID)
    fprintf(log);
else
    fprintf(fileID, log);
end % End if

end % End function
