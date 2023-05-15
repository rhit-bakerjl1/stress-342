%% Solving Hertzian Contact Stress
function [Resid] = func_vel_air_resid(v_init_mph, v_desired_mph, air_length)
    
    % Does v_desired exist?
    if (~exist("v_desired_mph", "var"))
        v_desired_mph   = 0;
    end

    if (~exist("air_length", 'var'))
        air_length  = 2*12;     % in
    end

    v_init  = v_init_mph*17.6; % in/s
    v_desired   = v_desired_mph*17.6;  % in/s

    % Plot?
    plotting    = 0;

    % Constants     (1 = baseball, 2 = ice)
    R   = 2.9/2;          % in
    m   = 0.3203;       % lbm, [0.3125, 0.3281]
    rho = 4.425593e-5;  % lbm/in^3
    Cd  = 0.3;          % unitless
    
    % Finally doing the baseball one
    d0      = 0;    % in
    v_vec0      = [d0, v_init];
    
    % Time
    t0  = 0;        % s
    % tf  = 0.0004;   % s
    tf  = 0.1;   % s
    N_t     = 3000*abs(log10(v_desired_mph))^2/abs(log10(1.5))^2;
    delT    = (tf-t0)/(N_t-1);
    
    % Solution
    [v, t, stopped] = func_rk4_modif_air(@(v)func_dvdt(v,m,rho,R,Cd), t0, v_vec0, delT, air_length);
        
    % Analysis
    dddt_vec    = v(:,2);
    
    % Chopping
    if (plotting)
        d_vec   = v(:,1);
        t       = t(dddt_vec > 0);
        d_vec   = d_vec(dddt_vec > 0);
    end
    dddt_vec    = dddt_vec(dddt_vec > 0);
    
    % Residual
    endVel  = dddt_vec(end);
    if (stopped)
        endVel  = 0;
    end
    Resid   = endVel - v_desired_mph*17.6;

    if (stopped)
        Resid   = -1000/v_init_mph - v_desired_mph*17.6;
    end

    % Plotting
    if (plotting)
        figure(1);
        clf;
        subplot(3,1,1);
        plot(t, d_vec);
        xlabel("Time (s)");
        ylabel("Displacement (in)");

        subplot(3,1,2);
        plot(t, dddt_vec);
        xlabel("Time (s)");
        ylabel("Velocity (in/s)");

        d2ddt2_vec  = (dddt_vec(2:end) - dddt_vec(1:end-1))/delT;
        subplot(3,1,3);
        plot(t(2:end),d2ddt2_vec);
        xlabel("Time (s)");
        ylabel("Acceleration (in/s^2)");                                                    
    end
end

%% Helpful Functions
function [dvdt] = func_dvdt(v, m, rho, R, Cd)
    % Constants
    A   = pi*R^2;   % in^2

    % Variables
    vel = v(2);     % in/s

    % Use the Force Luke
    F   = -rho*A*Cd*vel^2/2/m;   % in/s^2

    % Check if Stopped
    if (vel <= 0)
        vel  = 0;
        dvdt(2) = 0;
    end

    % Fill in dvdt
    dvdt(1) = vel;  % in/s
    dvdt(2) = F;    % in/s^2

end