% Computes the solution to a DE using RK4
% odefun: f(y,t) that computes y' in row vector form
% tspan: [t0, tf]
% y0: y @ t=t0. Row Vector
% delT: desired timestep
function [y, t] = func_rk4(odefun, tspan, y0, delT)
    % Time stuff
    t0  = tspan(1);
    tf  = tspan(2);
    t   = (t0:delT:tf)';

    % Initialize y
    y   = zeros(length(t), length(y0));
    y(1,:)  = y0;

    % Loop through to find y
    for i = 1:length(t)-1
        t_half  = (t(i+1)-t(i))/2;
        % First stage
        dy0     = odefun(y(i,:), t(i));
        y1      = y(i,:) + delT/2*dy0;
        % Second Stage
        dy1     = odefun(y1, t_half);
        y2      = y(i,:) + delT/2*dy1;
        % Third Stage
        dy2     = odefun(y2, t_half);
        y3      = y(i,:) + delT*dy2;
        % Fourth Stage
        dy3     = odefun(y3, t_half);
        y(i+1,:)    = y(i,:) + delT/6*(dy0 + 2*dy1 + 2*dy2 + dy3);
    end

end