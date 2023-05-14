% Computes the solution to a DE using a modified RK4
% Stops when we have gone through a 4-inch ice block
% odefun: f(y,t) that computes y' in row vector form. Also gives a tag
% which is 1 when we have hit maximum pressure and 0 when we have not.
% tspan: [t0, tf]
% y0: y @ t=t0. Row Vector
% delT: desired timestep
function [y, t, stopped] = func_rk4_modif_air(odefun, t0, y0, delT)
    % Time stuff
    N_alloc     = 100000;
    t   = (t0:delT:delT*N_alloc)';
    tag = 0;
    stopped = 0;

    % Initialize y
    y   = zeros(length(t), length(y0));
    y(1,:)  = y0;
    thick   = 24;

    % Loop through to find y
    i   = 1;
    while (y(i,1) < thick && stopped ~= 1)
        % First stage
        [dy0]     = odefun(y(i,:));
        y1      = y(i,:) + delT/2*dy0;
        % Second Stage
        [dy1]     = odefun(y1);
        y2      = y(i,:) + delT/2*dy1;
        % Third Stage
        [dy2]     = odefun(y2);
        y3      = y(i,:) + delT*dy2;
        % Fourth Stage
        [dy3]     = odefun(y3);
        delY    = delT/6*(dy0 + 2*dy1 + 2*dy2 + dy3);
        y(i+1,:)    = y(i,:) + delY;
        if (delY(2) == 0)
            stopped     = 1;
        end
        i = i + 1;
    end

end
