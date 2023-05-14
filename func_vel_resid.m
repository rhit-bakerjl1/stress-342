%% Solving Hertzian Contact Stress
function [Resid] = func_vel_resid(v_init_mph, v_desired_mph, crk_len)
    
    % Does v_desired exist?
    if (~exist("v_desired_mph", "var"))
        v_desired_mph   = 0;
    end

    % Plot?
    plotting    = 0;

    % Constants     (1 = baseball, 2 = ice)
    R   = 2.9;          % in
    E1  = 10*10^3;      % psi
    E2  = 1320*10^3;    % psi
    nu2 = 0.33;         % unitless
    m1  = 0.3203;       % lbm, [0.3125, 0.3281]
    nu1 = 0.3806;   % unitless, [0.3388,0.4225]
    
    % Meshing Decisions
    Nr  = 300;
    
    % Beginning maths
    E_star  = 1./((1-nu1.^2)/E1 + (1-nu2.^2)/E2);   % psi
        
    % Finally doing the baseball one
    d0      = 0;    % in
    % v_init_mph  = 200;  % mph
    v_init  = v_init_mph*17.6; % in/s
    v_vec0      = [d0, v_init];
    
    % Time
    t0  = 0;        % s
    % tf  = 0.0004;   % s
    tf  = 0.001;   % s
    N_t     = 3000*abs(log10(v_desired_mph))/abs(log10(1.5));
    delT    = (tf-t0)/(N_t-1);
    
    % Solution
    if (exist('crk_len', 'var'))
        [v, t, stopped] = func_rk4_modif(@(v)func_dvdt(v,m1,E_star(1),R, crk_len), t0, v_vec0, delT);
    else
        [v, t, stopped] = func_rk4_modif(@(v)func_dvdt(v,m1,E_star(1),R), t0, v_vec0, delT);
    end
    
    % Analysis
    d_vec   = v(:,1);
    dddt_vec    = v(:,2);
    
    % Chopping
    t       = t(dddt_vec > 0);
    d_vec   = d_vec(dddt_vec > 0);
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
function [dvdt, tag] = func_dvdt(v, m, E_star, R, crk_len)
    % Conversions
    lbf2base    = 32.17;    % ft*lbm/lbf/s^2
    ft2in       = 12;       % in/ft
    tag         = 0;

    % Variables
    d   = v(1);     % in
    vel = v(2);     % in/s

    % Nr
    Nr  = 100;

    % Finding force information
    [F, P_profile, ~]   = func_force(d, E_star, R, Nr);     % lbf, psi
    sigmaMax    = max(P_profile);  % psi
    sigmaBreak  = 870.2;   % psi
    if (exist("crk_len", 'var'))
        Kic     = 182.149;  % psi*sqrt(in)
        sigmaBreak2     = Kic/sqrt(pi*crk_len);
        sigmaBreak  = min(sigmaBreak, sigmaBreak2);
    end

    % Differential
    if (vel <= 0)
        vel  = 0;
        dvdt(2) = 0;
    else
        if (sigmaMax > sigmaBreak)
            F   = R^2*(pi*sigmaBreak)^3/(6*E_star^2);   % lbf
            tag = 1;
        end
        dvdt(2) = -F/m*lbf2base*ft2in;     % in/s^2
    end
    dvdt(1) = vel;  % in/s

end

% For timestepping oh boy yeah less go 
function [F_vec, P_mat, r_mat] = func_force(d_vec, E_star, R, Nr)
    % d_vec is a N_t by 1 vector: in
    % F will be a N_t by 1 vector: lbf
    % P_mat will be a N_t by Nr vector: psi
    % r_mat will be a N_t by Nr vector: in

    % Find the force
    N_t     = length(d_vec);
    F_vec   = 4/3*E_star*sqrt(R)*sqrt(d_vec.^3);    % lbf
    F_mat   = repmat(F_vec, 1, Nr);     % N_t by Nr; lbf
    a_vec   = sqrt(R*d_vec);            % in
    a_mat   = repmat(a_vec, 1, Nr);     % N_t by Nr; in

    % Filling in r_mat
    r_mat   = zeros(N_t, Nr);       % in
    for i = 1:N_t
        r_mat(i,:)  = linspace(0,a_vec(i),Nr);  % in
    end

    % Use the force
    p0_mat  = 3*F_mat./(2*pi*a_mat.^2);     % N_t by Nr; psi
    P_mat   = sqrt(1-r_mat.^2./a_mat.^2).*p0_mat;   % psi
end