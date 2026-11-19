function p = cp_params()
%CP_PARAMS  Default parameters for the cardioplegia heat-exchanger model.
%
%   p = cp_params() returns a struct with every geometric, physical and
%   numerical parameter used by CP_SOLVE.  SI units throughout, except
%   temperatures which are in degC (only differences enter the physics).
%
%   >>> REPLACE THE GEOMETRY / FLOW VALUES BELOW WITH THE REAL DEVICE DATA <<<
%   The defaults describe a generic tube-bundle cardioplegia heat exchanger
%   (blood inside N parallel tubes, water on the shell side) of the size
%   typically used clinically (~0.03 m^2 area, 300 mL/min blood).

%% ---------------- Geometry (per tube, all tubes identical) -------------
p.Ntube   = 40;          % number of parallel blood tubes            [-]
p.R0      = 0.75e-3;     % clean inner radius of a blood tube        [m]
p.L       = 0.15;        % heat-exchange length                       [m]
p.t_wall  = 0.15e-3;     % tube wall thickness                        [m]
p.k_wall  = 16;          % wall conductivity (stainless steel 316)    [W/m/K]

%% ---------------- Blood flow -------------------------------------------
p.Q_total = 300e-6/60;   % total blood (cardioplegia) flow, 300 mL/min [m^3/s]

%% ---------------- Blood thermophysical properties ----------------------
p.rho = 1060;            % density                                    [kg/m^3]
p.cp  = 3600;            % specific heat                              [J/kg/K]
p.k   = 0.50;            % thermal conductivity                       [W/m/K]

%% ---------------- Rheology ---------------------------------------------
% 'newtonian' | 'carreau-yasuda' | 'powerlaw'
p.rheology = 'newtonian';
p.mu_N     = 3.5e-3;     % Newtonian viscosity (Hct ~ 40%, 37 C)      [Pa s]
% Carreau-Yasuda parameters for whole blood (Boyd et al. 2007,
% after Cho & Kensey 1991):  mu = mu_inf + (mu0-mu_inf)[1+(lam*g)^a]^((n-1)/a)
p.CY = struct('mu0',0.056,'muinf',0.00345,'lambda',1.902,'n',0.22,'a',1.25);
% Power-law (used only for validation / comparison):  mu = K g^(n-1)
p.PL = struct('K',0.017,'n',0.708);

%% ---------------- Temperature-dependent viscosity (case i) -------------
% mu(g,T) = aT(T) * mu(g),   aT = exp(-beta*(T - Tref))
% beta ~ 0.02 1/K  (blood viscosity rises ~2 % per degC of cooling)
p.Tdep    = false;
p.beta_mu = 0.021;       % [1/K]
p.T_ref   = 37;          % [degC]

%% ---------------- Viscous dissipation (case h) -------------------------
p.viscDiss = false;

%% ---------------- Thermal boundary condition ---------------------------
% 'T'       : constant wall temperature (case a)
% 'q'       : constant wall heat flux (case b), q>0 = heat leaves blood
% 'coupled' : wall temperature from coupled coolant energy balance (case d)
p.bc     = 'T';
p.T_in   = 42;           % blood inlet temperature                    [degC]
p.T_wall = 32;           % wall temperature for bc='T'                [degC]
p.q_wall = 3000;         % wall heat flux for bc='q' (+ = cooling)    [W/m^2]

%% ---------------- Coolant (water) circuit, bc = 'coupled' --------------
p.flow_arrangement = 'counter';   % 'counter' | 'co'
p.Qc_total = 2/60/1000;  % total water flow, 2 L/min                  [m^3/s]
p.Tc_in    = 28;         % water inlet temperature                    [degC]
p.rho_c    = 998;        % water density                              [kg/m^3]
p.cp_c     = 4182;       % water specific heat                        [J/kg/K]
p.h_c      = 1500;       % water-side heat-transfer coefficient       [W/m^2/K]

%% ---------------- Transient model only (cp_solve_transient) ------------
p.rho_w   = 8000;        % wall density (stainless steel)             [kg/m^3]
p.c_w     = 500;         % wall specific heat                         [J/kg/K]
p.wallCap = true;        % include wall heat capacity
p.Ac      = 3e-6;        % shell-side water flow area per tube        [m^2]

%% ---------------- Wall deposit / clot (case g) -------------------------
% Gaussian deposit:  R(x) = R0 * [1 - delta*exp(-((x-xc)/w)^2)]
p.clot.on    = false;
p.clot.delta = 0.40;     % max fractional radius reduction            [-]
p.clot.xc    = 0.50;     % centre, as fraction of L                   [-]
p.clot.w     = 0.010;    % half-width                                 [m]
p.clot.k     = 0.45;     % conductivity of thrombus layer             [W/m/K]

%% ---------------- Safety limits ----------------------------------------
p.T_safe      = 40;      % max cup-mixing blood temp (hemolysis)      [degC]
p.T_target    = 37;      % target outlet temperature                  [degC]
p.T_protein   = 42;      % local protein-denaturation limit           [degC]
p.p_out_mmHg  = 80;      % pressure downstream of the HX (delivery)   [mmHg]
p.p_max_mmHg  = 150;     % allowed max line pressure                  [mmHg]
p.tau_max     = 150;     % wall-shear-stress limit                    [Pa]
% Giersiepen et al. (1990) hemolysis power law  HI[%] = C t^alpha tau^beta
p.hem = struct('C',3.62e-5,'alpha',0.785,'beta',2.416);

%% ---------------- Numerics ---------------------------------------------
p.Nr         = 101;      % radial nodes (clustered at the wall)
p.Nx         = 500;      % axial steps   (clustered at the inlet)
p.x_cluster  = 2;        % x = L*s^x_cluster
p.picard_tol = 1e-7;     % [K]  Picard tolerance (T-dependent viscosity)
p.picard_max = 25;
p.shoot_tol  = 1e-6;     % [K]  counter-current shooting tolerance
end
