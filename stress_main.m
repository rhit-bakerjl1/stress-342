close;
clear;
clc;

% Options
onlyOne     = 1;

% Main function
if (onlyOne)
    options     = ones(1,4);
    % Initial Velocity
    v0  = 200;  % mph
    func_disp_press(v0, options);
else
    options     = zeros(1,4);
    % Initial Velocity
    v0_min  = 0.1;   % mph
    v0_max  = 1000;  % mph
    N_v     = 400;  
    v0_vec  = linspace(v0_min, v0_max, N_v);

    % Iterating through for each velocity
    max_d_vec   = zeros(size(v0_vec));
    max_P_vec   = zeros(size(v0_vec));
    max_d_exact = zeros(size(v0_vec));
    for i = 1:N_v
        [d_vec, P_mat, max_d_exact(i)]  = func_disp_press(v0_vec(i), options);
        max_d_vec(i)    = max(d_vec);
        max_P_vec(i)    = max(P_mat, [], 'all');
    end

    % Plotting
    figure(1);
    clf;
    plot(v0_vec, max_d_vec, "--");
    hold on;
    plot(v0_vec, max_d_exact);
    xlabel("Initial Velocity (mph)");
    ylabel("Maximum deflection (in)");
    legend("Approximate", "Exact Model");

    figure(2);
    clf;
    plot(v0_vec, max_P_vec);
    xlabel("Initial Velocity (mph)");
    ylabel("Maximum Pressure (ksi)");

end

