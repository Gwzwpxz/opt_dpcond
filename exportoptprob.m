function [sprob] = exportoptprob(prob, path, fname)
% Export an optimal preconditioning problem to sedumi file and write it to
% local SDPA file. The SDPA writer comes from CSDP.

writetofile = true;
if nargin < 3
    writetofile = false;
end % End if

ptype = upper(prob.ptype);

if ptype == 'R'
    [A, b, c, K] = getsedumi(prob.M, ptype);
elseif ptype == 'L'
    [A, b, c, K] = getsedumi(prob.X, ptype);
else
    error("Only 'R' and 'L' problems can be exported to single SDP");
end % End if

sprob.A = A;
sprob.b = b;
sprob.c = c;
sprob.K = K;

if writetofile
    param.printlevel = 0;
    sdpapath = fullfile(path, string(fname) + "-" + ptype + ".dat-s");
    writesdpa(sdpapath, A, b, c, K, param);
    fprintf("Exporting SDP to %s \n", sdpapath);
end % End if

end % End function