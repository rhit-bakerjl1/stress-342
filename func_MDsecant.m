function [soln, er_est] = func_MDsecant(R, xi, tol, maxIter, toggle)
    count = 0;
    dx = tol * 50;
    x = xi;
    f = R(xi);
    x2      = 1.01*x;
    delx    = 0.01*x;

    f2  = R(x2);
    df  = (f2-f)/delx;
    if (toggle)
        fprintf('\n Count       xi           R(xi)         dR(xi)        dx \n');
    end
    while(abs(dx) > tol && count < maxIter)
        count   = count + 1;
        dx      = -f/df;
        if (toggle)
            fprintf('   %3i     % 10.3e   % 10.3e    % 10.3e    % 10.3e\n', count, x, f, df, dx);
        end
        x       = x + dx;
        f       = R(x);
        x2      = 1.01*x;
        delx    = 0.01*x;

        f2  = R(x2);
        df  = (f2-f)/delx;
    end
    soln = x;
    er_est = abs(f/df);
end