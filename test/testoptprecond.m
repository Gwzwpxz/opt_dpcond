function [D] = testoptprecond(data, type, solver, pcond, testsolve, fileID)

if nargin == 5
    fileID = [];
end % End if

if nargin == 4
    testsolve = false;
    fileID = [];
end % End if

DSDP_RETCODE_OK = 2;
sdptime = 0.0;

tpcg = -1.0;
tdirect = -1.0;
itpcg = -1;

if pcond == "L"
    M = data.M;
    cond_before = cond(full(M), 2);
    if solver == "cvx"
        D = getcvxdiag(M, pcond);
    else
        % Apply command-line tool of HDSDP to get solution
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
    
    PM = sqrt(D) \ (sqrt(D) \ full(data.M))';
    cond_after = cond(PM);
    reduce_cond = cond_before - cond_after;
    
elseif pcond == "R"
    
    % GET RHS pre-conditioning
    X = data.X; XTX = X' * X;
    cond_before = cond(full(XTX));
    if cond_before > 1e+08
        return;
    else
        if solver == "cvx"
            D = getcvxdiag(X, pcond);
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
            catch
                return;
            end % End try
        end % End if
        
        XDX = X' * D * X;
        cond_after = cond(full(XDX));
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
        DXTXD = sqrt(d1) \ (sqrt(d1) \ XTX)';
        cond_after = cond(DXTXD);
    catch
        return;
    end % End try
    reduce_cond = cond_before - cond_after;
    
end % End if

% Numerical issue happens
if reduce_cond < 0
    return;
end % End if

% Check how good optimal pre-conditioner is 
rdc = reduce_cond / cond_before;
class = 3;
if rdc >= 0.8
    class = 1;
elseif rdc >= 0.5
    class = 2;
end % End if

if class <= 2
    
    if pcond == "L"
        % Compare with diagonal preconditioner
        Ddiag = diag(diag(data.M));
        MDiag = sqrt(Ddiag) \ (sqrt(Ddiag) \ full(data.M))';
        cond_jacob = cond(MDiag);
        
        if testsolve
            [tpcg, tdirect, ~, itpcg] = getslvtime(XTX, XDX, 12800);
        end % End if
        
        log = sprintf("%d(%s) %30s %6d | %6.3e %6.3e (J: %6.3e) %f T: %f || %3.3f(%d) %3.3f ||| %d %d %d \n",...
            class, pcond, data.name, size(D, 1), cond_before,...
            cond_after, cond_jacob, reduce_cond / cond_before,...
            sdptime, tpcg, itpcg, tdirect);
        
    elseif pcond == "R"
        
        if testsolve
            [tpcg, tdirect, ~, itpcg] = getslvtime(XTX, XDX, 12800);
        end % End if 
        log = sprintf("%d(%s) %30s %6d | %6.3e %6.3e %f T: %f || %3.3f(%d) %3.3f \n", class, pcond, data.name,...
            size(D, 1), cond_before, cond_after, reduce_cond / cond_before, sdptime, tpcg, itpcg, tdirect);
    else
        
        if testsolve
            [tpcg, tdirect, ~, itpcg] = getslvtime(XTX, DXTXD, 128000);
        end % End if
        
        [L, R] = ruizscale(data.X, 100);
        Xcond = diag(L) * data.X * diag(R);
        cond_ruiz = cond(full(Xcond' * Xcond));
        
        log = sprintf("%d(%s) %30s %6d | %6.3e %6.3e %6.3e T: %f || %3.3f(%d) %3.3f\n",...
            class, pcond, data.name, size(X, 1), cond_before,...
            cond_after, cond_ruiz, reduce_cond / cond_before,...
            tpcg, itpcg, tdirect);
        
    end % End if
    
end % End if

if isempty(fileID)
    fprintf(log);
else
    fprintf(fileID, log);
end % End if

end % End function
