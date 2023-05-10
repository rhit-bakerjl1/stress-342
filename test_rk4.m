%% Testing RK4 with 1D Solution
close;
clear;
clc;

% Constants
t0  = 0;
tf  = 3;
delT    = 1;
y0  = 3;
alpha   = 2;

% Finding Approximate Solution
[v,t] = func_rk4(@(y,t)func_dydt(y,t,alpha), [t0, tf], y0, delT);

% Finding Exact Solution
y_exact     = y0*exp(alpha*t);

% Error
err_vec     = y_exact - v;

% Plotting Result
figure(1);
clf;
plot(t, v, 'o');
hold on;
plot(t, y_exact);
xlabel("Time (s)");
ylabel("y");
legend("RK4 Solution", "Exact Solution");

% Plotting Error
figure(2);
clf;
plot(t, err_vec, "o");
xlabel("Time (s)");
ylabel("Error");

% Displaying Results
disp("Maximum Error is: " + max(err_vec));

%% Testing RK4 with 2D Solution
% Constants
t0  = 0;
tf  = 3;
delT    = 0.2;
y0  = 0;
beta    = 3;
dydt0   = 2;
A       = dydt0/beta;

% Finding Approximate Solution
[v,t] = func_rk4(@(y,t)func_dvdt(y,t,beta), [t0, tf], [y0; dydt0], delT);
y_aprx  = v(:,1);
dydt_aprx  = v(:,2);

% Finding Exact Solution
y_exact     = A*sin(beta*t);
dydt_exact  = A*beta*cos(beta*t);

% Error
err_vec     = y_exact - y_aprx;

% Plotting Result
figure(3);
clf;
plot(t, y_aprx, 'o');
hold on;
plot(t, y_exact);
xlabel("Time (s)");
ylabel("y");
legend("RK4 Solution", "Exact Solution");

figure(4);
clf;
plot(t, dydt_aprx, 'o');
hold on;
plot(t, dydt_exact);
xlabel("Time (s)");
ylabel("dy/dt");
legend("RK4 Solution", "Exact Solution");

% Plotting Error
figure(5);
clf;
plot(t, err_vec, "o");
xlabel("Time (s)");
ylabel("Error");

% Displaying Results
disp("Maximum Error is: " + max(err_vec));

%% Functions
function [dydt] = func_dydt(y, t, alpha)
    dydt    = alpha*y;
end

function [dvdt] = func_dvdt(v, t, beta)
    v1  = v(1);
    v2  = v(2);

    dvdt(1)     = v2;
    dvdt(2)     = -beta^2*v1;
end