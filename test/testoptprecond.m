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
    
elseif pcond == "M"
    % Mixed pre-conditioner
    M = full(data.M);
    
%     fprintf("%30s Sp: %6.3e \n", data.name, cond_spinv);
%     return;
    
    cond_before = cond(full(M), 2);
    
    % Two base pre-conditioners
    Ddiag = diag(diag(M));
    [Druiz, Dtmp] = ruizscale(M, 100); %#ok
    MDiag = sqrt(Ddiag) \ (sqrt(Ddiag) \ full(data.M))';
    Mruiz = Druiz * M * Druiz;
%     Dspinv = diag(diag(M) ./ sum(M.^2, 2));
%     Mspinv = sqrt(Dspinv) * M * sqrt(Dspinv);
%     cond_spinv = cond(Mspinv);
    cond_jacob = cond(MDiag);
    cond_ruiz = cond(Mruiz);
    
    Druiz = diag(1./diag(Druiz.^2));
%     Dspinv = diag(1./diag(Dspinv));
    
    try
        [D, cvx_time] = getoptcombprecond(M, Ddiag, Druiz, eye(size(M, 1)));
        if min(diag(D)) < 0
            error("!");
        end % End if 
    catch
        return;
    end % End try
    
    Mcomb = sqrt(D) \ (sqrt(D) \ M)';
    
    try 
        cond_comb = cond(Mcomb);
    catch
        return;
    end % End try
    
    bench = min([cond_jacob, cond_ruiz, cond_before]);
    
    factor = bench / cond_comb;
    
    class = 3;
    if factor >= 2.0
        class = 1;
    elseif factor >= 1.2
        class = 2;
    end % End if
        
    % Numerical issue happens
    if factor <= 1.0
        return;
    end % End if
    
    % Test conjugate gradient
    if testsolve && class <= 2
        
        tol = 1e-04;
        maxit = size(Mcomb, 1) * 10;
        rng(maxit);
        
        rhs = randn(size(Mcomb, 1), 1);
        
        Mcomb = sparse(Mcomb) / trace(Mcomb);
        Mruiz = sparse(Mruiz) / trace(Mruiz);
        MDiag = sparse(MDiag) / trace(MDiag);
        
        [xcomb, flagcomb, rescomb, itcomb] = pcg(Mcomb, sqrt(D) \ rhs, tol, maxit); %#ok
        [xruiz, flagruiz, resruiz, itruiz] = pcg(Mruiz, sqrt(Druiz) \ rhs, tol, maxit); %#ok
        [xjacob, flagjacob, resjacob, itjacob] = pcg(MDiag, sqrt(Ddiag) \ rhs, tol, maxit); %#ok
        % [xsp, flagsp, ressp, itsp] = pcg(Mspinv, sqrt(Dspinv) \ rhs, tol, maxit); %#ok
        [xorig, flagorig, resorig, itorig] = pcg(M, rhs, tol, maxit); %#ok
        
        if flagcomb ~= 0 
            itcomb = maxit;
        end % End if 
        
        if flagruiz ~= 0 
            itruiz = maxit;
        end % End if 
        
        if flagjacob ~= 0 
            itjacob = maxit;
        end % End if 
        
%         if flagsp ~= 0 
%             itsp = maxit;
%         end % End if 
        
        if flagorig ~= 0 
            itorig = maxit;
        end % End if 
        
        if itcomb == min([itcomb, itruiz, itjacob, itorig])
            status = "*";
        else
            status = "x";
        end % End if
        
        fprintf("%d(%s) %30s          O: %4d | J: %4d | R: %4d | C: %4d | %s \n",...
            class, "M", data.name, itorig, itjacob, itruiz, itcomb, status);
        
    end % End if
    
    log = sprintf("%d(%s) %30s %6d | O: %6.3e J: %6.3e R: %6.3e C: %6.3e || t: %f facmin: %f facimp: %f\n", class, "M",...
        data.name, size(D, 1), cond_before, cond_jacob, cond_ruiz, cond_comb, cvx_time, factor - 1, cond_before / cond_comb - 1);
    
    if isempty(fileID)
        fprintf(log);
    else
        fprintf(fileID, log);
    end % End if
    
    return;
    
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
factor = cond_before / cond_after;
class = 3;
if factor >= 5.0
    class = 1;
elseif factor >= 2.0 
    class = 2;
end % End if

if class <= 3
    
    if pcond == "L"
        % Compare with diagonal preconditioner
        Ddiag = diag(diag(data.M));
        MDiag = sqrt(Ddiag) \ (sqrt(Ddiag) \ full(data.M))';
        cond_jacob = cond(MDiag);
        
        if testsolve
            [tpcg, tdirect, ~, itpcg] = getslvtime(XTX, XDX, 12800);
        end % End if
        
        log = sprintf("%d(%s) %30s %6d | %6.3e %6.3e (J: %6.3e) %f T: %f || %3.3f(%d) %3.3f \n",...
            class, pcond, data.name, size(D, 1), cond_before,...
            cond_after, cond_jacob, factor,...
            sdptime, tpcg, itpcg, tdirect);
        
    elseif pcond == "R"
        
        if testsolve
            [tpcg, tdirect, ~, itpcg] = getslvtime(XTX, XDX, 12800);
        end % End if 
        log = sprintf("%d(%s) %30s %6d | %6.3e %6.3e %f T: %f || %3.3f(%d) %3.3f \n",...
            class, pcond, data.name, size(D, 1), cond_before,...
            cond_after, factor,...
            sdptime, tpcg, itpcg, tdirect);
    else
        
        if testsolve
            [tpcg, tdirect, ~, itpcg] = getslvtime(XTX, DXTXD, 128000);
        end % End if
        
        [L, R] = ruizscale(data.X, 100);
        Xcond = diag(L) * data.X * diag(R);
        cond_ruiz = cond(full(Xcond' * Xcond));
        
        log = sprintf("%d(%s) %30s %6d | %6.3e %6.3e %6.3e T: %f || %3.3f(%d) %3.3f\n",...
            class, pcond, data.name, size(X, 1), cond_before,...
            cond_after, cond_ruiz, factor,...
            tpcg, itpcg, tdirect);
        
    end % End if
    
end % End if

if isempty(fileID)
    fprintf(log);
else
    fprintf(fileID, log);
end % End if

end % End function
