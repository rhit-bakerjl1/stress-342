%% Solving Hertzian Contact Stress
close;
clear;
clc;

% Options
test_nu1    = 1;

% Constants     (1 = baseball, 2 = ice)
R   = 2.9;          % in
E1  = 10*10^3;      % psi
E2  = 1320*10^3;    % psi
nu2 = 0.33;         % unitless
m1  = 0.3203;       % lbm, [0.3125, 0.3281]
if (test_nu1)
    nu1_min     = 0;
    nu1_max     = 0.7;
    N_nu        = 200;
    nu1     = linspace(nu1_min,nu1_max,N_nu);  % unitless
else
    nu1     = 0.3806;   % unitless, [0.3388,0.4225]
end

% Meshing Decisions
Nr  = 300;
d_test  = 0.5;  % in

% Plotting number
figNum  = 1;

% Conversions
lbf2base     = 32.17;    % ft*lbm/lbf/s^2
ft2in        = 12;       % in/ft

% Beginning maths
E_star  = 1./((1-nu1.^2)/E1 + (1-nu2.^2)/E2);   % psi
alpha   = 4*E_star*sqrt(R)*lbf2base*ft2in/3/m1; % 1/(sqrt(in)*s^2)

% Does Poisson's Ratio Matter?
if (test_nu1)
    figure(figNum);
    figNum = figNum + 1;
    clf;
    % plot(nu1, E_star);
    plot(nu1, alpha);
    xlabel("Baseball Poisson's Ratio nu");
    % ylabel("E_star");
    ylabel("Alpha");
end

% Find Pressure and Stuff at given d
[F, P_vec, r_vec]   = func_force_E_star(d_test, E_star, R, Nr);

% Plotting Pressure distribution
figure(figNum);
figNum  = figNum + 1;
clf;
r_vec_plt   = [-flip(r_vec(2:end)); r_vec];
if (test_nu1)
    % Getting three different nu values to check
    P_plt_ind   = ceil(linspace(1, N_nu, 5));
    hold on;
    legendStrs  = strings(1,length(P_plt_ind));
    for i = 1:length(P_plt_ind)
        ind     = P_plt_ind(i);
        P_vec_plt   = [flip(P_vec(2:end,ind)); P_vec(:,ind)];
        plot(r_vec_plt, P_vec_plt);
        legendStrs(i)   = "nu = " + nu1(ind);
    end
    % Plotting pressure for these nu values
    legend(legendStrs);
else
    % Plot
    P_vec_plt   = [flip(P_vec(2:end)); P_vec];
    plot(r_vec_plt, P_vec_plt);
end
xlabel("Radius from center of pressure circle (in)");
ylabel("Pressure (psi)");

% Finally doing the baseball one
d0      = 0;    % in
v_init_mph  = 200;  % mph
v_init  = v_init_mph*17.6; % in/s
v_vec0      = [d0, v_init];

% Time
t0  = 0;        % s
tf  = 0.0004;   % s
N_t     = 300;
delT    = (tf-t0)/(N_t-1);

% Solution
[v, t] = func_rk4(@(v,t)func_dvdt(v,t,alpha(1)), [t0, tf], v_vec0, delT);

% Analysis
d_vec   = v(:,1);
dddt_vec    = v(:,2);

% Chopping
t       = t(dddt_vec > 0);
d_vec   = d_vec(dddt_vec > 0);
dddt_vec    = dddt_vec(dddt_vec > 0);

% Finding pressures
[F_vec, P_mat, r_mat] = func_force(d_vec, E_star(1), R, Nr);

% Pressure Movie!
figure(figNum);
figNum = figNum + 1;
clf;
% Min/max
minR    = -max(r_mat, [], 'all');
maxR    = -minR;
minP    = min(P_mat, [], 'all');
maxP    = max(P_mat, [], 'all');
% Iteration
for i = 2:length(t)
    % Figure
    figure(figNum-1);
    % Extending r_plot and P_plot
    r_plot  = [-flip(r_mat(i,2:end)), r_mat(i,:)];
    P_plot  = [flip(P_mat(i,2:end)), P_mat(i,:)];
    % Plotting
    plot(r_plot, P_plot);
    xlabel("Distance from center (in)");
    ylabel("Pressure (ksi)");
    xlim([minR*1.1, maxR*1.1]);
    ylim([minP*1.1, maxP*1.1]);
    title("Time = " + round(t(i),6) + " s");
    pause(1/48);
end

% Deflection Plotting
figure(figNum);
figNum = figNum + 1;
clf;
plot(t, d_vec);
xlabel("Time (s)");
ylabel("Deflection into material (in)");

%% Helpful Functions
function [dvdt] = func_dvdt(v, t, alpha)
    % Variables
    x1  = v(1);
    x2  = v(2);

    % Differential
    if (x2 <= 0)
        x2  = 0;
        dvdt(2) = 0;
    else
        dvdt(2) = -alpha*x1^(3/2);
    end
    dvdt(1) = x2;

end

% For testing if E_star is important or not
function [F_vec, P_vec, r_vec] = func_force_E_star(d, E_star, R, Nr)
    % Find the force
    F_vec   = 4/3*E_star*sqrt(R)*sqrt(d^3);
    a       = sqrt(R*d);

    % Filling in r_mat
    r_vec   = linspace(0,a,Nr)';

    % Use the force
    p0_vec  = 3*F_vec./(2*pi*a.^2);     % N_t by Nr
    P_vec   = sqrt(1-r_vec.^2/a^2)*p0_vec;
end

% For timestepping oh boy yeah less go 
function [F_vec, P_mat, r_mat] = func_force(d_vec, E_star, R, Nr)
    % d_vec is a N_t by 1 vector
    % F will be a N_t by 1 vector
    % P_mat will be a N_t by Nr vector
    % r_mat will be a N_t by Nr vector

    % Find the force
    N_t     = length(d_vec);
    F_vec   = 4/3*E_star*sqrt(R)*sqrt(d_vec.^3);    % lbf
    F_vec   = F_vec / 1000;             % kip
    F_mat   = repmat(F_vec, 1, Nr);     % N_t by Nr
    a_vec   = sqrt(R*d_vec);
    a_mat   = repmat(a_vec, 1, Nr);     % N_t by Nr

    % Filling in r_mat
    r_mat   = zeros(N_t, Nr);
    for i = 1:N_t
        r_mat(i,:)  = linspace(0,a_vec(i),Nr);
    end

    % Use the force
    p0_mat  = 3*F_mat./(2*pi*a_mat.^2);     % N_t by Nr
    P_mat   = sqrt(1-r_mat.^2./a_mat.^2).*p0_mat;
end
