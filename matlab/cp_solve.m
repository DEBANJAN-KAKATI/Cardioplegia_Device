function S = cp_solve(p)
%CP_SOLVE  Steady 2-D (axisymmetric r-x) laminar conjugate heat transfer of
%          blood in one tube of the cardioplegia heat exchanger.
%
%   Energy equation (lubrication approximation, axial conduction neglected
%   since Pe = Re*Pr >> 1):
%
%     rho cp ( u dT/dx + v dT/dr ) = (1/r) d/dr( k r dT/dr ) + mu(gd,T) gd^2
%
%   Solved in the mapped coordinate eta = r/R(x) (handles wall deposits),
%   where the convective operator becomes  u T_x - (G_x /(R^2 eta)) T_eta,
%   G = R^2 int_0^eta u eta' d eta' (stream function / 2pi).
%
%   Wall conditions (p.bc):
%     'T'       T(R0) = Tw                               (case a)
%     'q'       -k dT/dr = q_w                           (case b)
%     'coupled' -k dT/dr = U (T_wall - Tc(x)),
%               mc cpc dTc/dx = +/- 2 pi R U (T_wall - Tc)  (case d)
%               counter-current solved by shooting on Tc(x=0).
%   With a wall deposit (p.clot.on) the thrombus layer adds a conduction
%   resistance R ln(R0/R)/k_clot in series.
%
%   Discretisation: finite volume in eta (wall-clustered), fully implicit
%   marching in x (inlet-clustered), Picard iteration when viscosity
%   depends on temperature.
%
%   S = cp_solve(p) returns fields along x (Tm, T_wall, q_w, Nu, p, tau_w,
%   Tc, ...) and the full T(eta,x), u(eta,x), mu(eta,x) fields.

Nt = p.Ntube;
Qt = p.Q_total/Nt;                       % flow per tube

% ---------- grids
s    = linspace(0,1,p.Nr)';
eta  = sin(pi/2*s);                      % clustered at the wall eta = 1
x    = p.L*linspace(0,1,p.Nx+1).^p.x_cluster;
g.eta = eta;  g.x = x;
g.h    = diff(eta);
g.etaf = 0.5*(eta(1:end-1)+eta(2:end));
g.V    = ([g.etaf;1].^2 - [0;g.etaf].^2)/2;   % CV "volumes"  int eta d eta

mcp_c = p.rho_c*p.cp_c*p.Qc_total/Nt;    % coolant capacity rate per tube

if strcmpi(p.bc,'coupled') && strcmpi(p.flow_arrangement,'counter')
    % ---- shooting: find coolant outlet temp Tc(0) s.t. Tc(L) = Tc_in
    res = @(Tc0) march(p, g, Qt, mcp_c, Tc0);
    a  = p.Tc_in;          Sa = res(a);  fa = Sa.Tc(end) - p.Tc_in;
    b  = p.Tc_in + 0.5*(p.T_in - p.Tc_in)*0.2 + 0.1;
    Sb = res(b);           fb = Sb.Tc(end) - p.Tc_in;
    for it = 1:30
        if abs(fb) < p.shoot_tol, break; end
        c  = b - fb*(b-a)/(fb-fa);
        a = b; fa = fb;
        b = c; Sb = res(b); fb = Sb.Tc(end) - p.Tc_in;
    end
    S = Sb;  S.shoot_iter = it;  S.shoot_res = fb;
else
    S = march(p, g, Qt, mcp_c, p.Tc_in);
end
S.p = p;
end

% ======================================================================
function S = march(p, g, Qt, mcp_c, Tc0)
eta = g.eta; x = g.x; N = numel(eta); Nx = numel(x);
rc  = p.rho*p.cp;
coupled = strcmpi(p.bc,'coupled');
sgn = 1; if coupled && strcmpi(p.flow_arrangement,'counter'), sgn = -1; end
Ntot = N + coupled;

% central first-derivative weights on the non-uniform grid (interior)
hm = g.h(1:end-1); hp = g.h(2:end);
Dm = -hp./(hm.*(hm+hp));  D0 = (hp-hm)./(hm.*hp);  Dp = hm./(hp.*(hm+hp));
% one-sided 2nd-order derivative at the wall (eta=1)
h1 = g.h(end); h2 = g.h(end-1);
Dw = [h1/(h2*(h1+h2)), -(h1+h2)/(h1*h2), (2*h1+h2)/(h1*(h1+h2))];  % on N-2,N-1,N

% storage
T   = p.T_in*ones(N,1);
R   = cp_radius(0,p);
Vel = cp_velocity(R, Qt, eta, T, p);
Gold = R^2*cumtrapz(eta, Vel.u.*eta);
Tc  = Tc0;

S.x = x;  S.eta = eta;
S.R = zeros(1,Nx);  S.Tm = S.R; S.Twall = S.R; S.Tc = S.R; S.qw = S.R;
S.Nu = S.R; S.tau_w = S.R; S.dpdx = S.R; S.Tmax = S.R; S.Tmin = S.R;
S.Tfield = zeros(N,Nx); S.ufield = S.Tfield; S.mufield = S.Tfield;
S.diss = S.R;                                   % dissipation per length [W/m]
S.qw_cons = S.R;

S.R(1) = R; S.Tm(1) = p.T_in; S.Twall(1) = p.T_in; S.Tc(1) = Tc;
S.qw(1) = NaN; S.Nu(1) = NaN; S.tau_w(1) = Vel.tau_w; S.dpdx(1) = Vel.dpdx;
S.Tmax(1) = p.T_in; S.Tmin(1) = p.T_in;
S.Tfield(:,1) = T; S.ufield(:,1) = Vel.u; S.mufield(:,1) = Vel.mu;
S.diss(1) = 2*pi*R^2*trapz(eta, Vel.Phi.*eta)*p.viscDiss;
S.qw_cons(1) = 0;

hemI = zeros(N,1);  ab = p.hem.beta/p.hem.alpha;
uold = Vel.u;  tauold = Vel.tau;
Rcache = R; Vcache = Vel;
dxv = [0 diff(x)];

for i = 2:Nx
    dx = x(i)-x(i-1);
    R  = cp_radius(x(i),p);
    [type, U, Text] = wall_bc(p, R, Tc);
    Tn = T;
    for it = 1:p.picard_max
        if ~p.Tdep && R == Rcache          % profile depends on R only -> reuse
            Vel = Vcache;
        else
            Vel = cp_velocity(R, Qt, eta, Tn, p);
            Rcache = R; Vcache = Vel;
        end
        G   = R^2*cumtrapz(eta, Vel.u.*eta);
        W   = zeros(N,1);
        W(2:N-1) = (G(2:N-1)-Gold(2:N-1))/dx ./ (R^2*eta(2:N-1));
        u   = Vel.u;

        % ---- assemble  A*[T;Tc] = b
        Ii = zeros(8*N,1); Jj = Ii; Vv = Ii; nz = 0; b = zeros(Ntot,1);
        % axial convection
        cd = rc*g.V.*u/dx;
        idx = (1:N)';  [Ii,Jj,Vv,nz] = add(Ii,Jj,Vv,nz, idx, idx, cd);
        b(1:N) = cd.*T;
        % mapped radial convection  -rc V W T_eta  (interior nodes)
        j  = (2:N-1)';  cw = -rc*g.V(j).*W(j);
        [Ii,Jj,Vv,nz] = add(Ii,Jj,Vv,nz, j, j-1, cw.*Dm);
        [Ii,Jj,Vv,nz] = add(Ii,Jj,Vv,nz, j, j,   cw.*D0);
        [Ii,Jj,Vv,nz] = add(Ii,Jj,Vv,nz, j, j+1, cw.*Dp);
        % radial conduction (faces)
        a  = p.k*g.etaf./g.h/R^2;  f = (1:N-1)';
        [Ii,Jj,Vv,nz] = add(Ii,Jj,Vv,nz, f,   f,    a);
        [Ii,Jj,Vv,nz] = add(Ii,Jj,Vv,nz, f,   f+1, -a);
        [Ii,Jj,Vv,nz] = add(Ii,Jj,Vv,nz, f+1, f+1,  a);
        [Ii,Jj,Vv,nz] = add(Ii,Jj,Vv,nz, f+1, f,   -a);
        % viscous dissipation
        if p.viscDiss, b(1:N) = b(1:N) + g.V.*Vel.Phi; end
        % wall
        switch type
            case 'neumann'
                b(N) = b(N) - U/R;                 % here U holds q_w
            case 'robin'
                [Ii,Jj,Vv,nz] = add(Ii,Jj,Vv,nz, N, N, U/R);
                b(N) = b(N) + U/R*Text;
            case 'robin_coupled'
                K = 2*pi*R*U/mcp_c;
                [Ii,Jj,Vv,nz] = add(Ii,Jj,Vv,nz, N,   N,    U/R);
                [Ii,Jj,Vv,nz] = add(Ii,Jj,Vv,nz, N,   N+1, -U/R);
                [Ii,Jj,Vv,nz] = add(Ii,Jj,Vv,nz, N+1, N+1, 1/dx + sgn*K);
                [Ii,Jj,Vv,nz] = add(Ii,Jj,Vv,nz, N+1, N,   -sgn*K);
                b(N+1) = Tc/dx;
        end
        A = sparse(Ii(1:nz),Jj(1:nz),Vv(1:nz),Ntot,Ntot);
        if strcmp(type,'dirichlet')
            A(N,:) = 0; A(N,N) = 1; b(N) = Text;  %#ok<SPRIX>
        end
        sol  = A\b;
        Tnew = sol(1:N);
        if ~p.Tdep, break; end
        if max(abs(Tnew-Tn)) < p.picard_tol, break; end
        Tn = Tnew;
    end
    T = Tnew;
    if coupled, Tc = sol(N+1); end

    % ---- post-processing at x(i)
    Tm = trapz(eta,u.*T.*eta)/trapz(eta,u.*eta);
    % conservative wall flux from the wall-cell balance (exact discrete energy)
    Fm  = p.k*g.etaf(end)*(T(N)-T(N-1))/g.h(end);
    Fw  = Fm - R^2*g.V(N)*Vel.Phi(N)*p.viscDiss;      % k*eta*dT/deta at wall
    qcons = -Fw/R;
    switch type
        case 'dirichlet', qw = -(p.k/R)*(Dw*T(N-2:N));   % 2nd-order gradient
        case 'neumann',   qw = U;
        case 'robin',     qw = U*(T(N)-Text);
        case 'robin_coupled', qw = U*(T(N)-Tc);
    end
    S.R(i) = R; S.Tm(i) = Tm; S.Twall(i) = T(N); S.Tc(i) = Tc; S.qw(i) = qw;
    S.qw_cons(i) = qcons;
    S.Nu(i) = 2*R*qw/(p.k*(Tm-T(N)));        % >0 for heat leaving or entering
    S.tau_w(i) = Vel.tau_w; S.dpdx(i) = Vel.dpdx;
    S.Tmax(i) = max(T); S.Tmin(i) = min(T);
    S.Tfield(:,i) = T; S.ufield(:,i) = u; S.mufield(:,i) = Vel.mu;
    S.diss(i) = 2*pi*R^2*trapz(eta, Vel.Phi.*eta)*p.viscDiss;

    % ---- hemolysis (Garon-Farinas linearised Giersiepen law along eta=const)
    jj = 1:N-1;
    hemI(jj) = hemI(jj) + 0.5*( tauold(jj).^ab./uold(jj) + Vel.tau(jj).^ab./u(jj) )*dx;
    uold = u; tauold = Vel.tau;
    Gold = G;
end

% ---------- integral quantities
S.p_rel = -cumtrapz(x, S.dpdx);                         % p(0)-p(x)
S.dp    = S.p_rel(end);                                  % pressure drop [Pa]
mmHg    = 133.322;
S.pres  = p.p_out_mmHg*mmHg + S.dp - S.p_rel;             % gauge pressure p(x) [Pa]
S.p_in_mmHg = S.pres(1)/mmHg;
S.dp_mmHg   = S.dp/mmHg;
HI_eta  = p.hem.C*hemI.^p.hem.alpha;                     % % per streamline
w       = S.ufield(:,end).*eta;
S.HI    = trapz(eta, HI_eta.*w)/trapz(eta, w);           % flow-averaged [%]
S.HI_eta = HI_eta;
% global energy balance per tube
mcp_b = p.rho*p.cp*Qt;
E_b   = mcp_b*(S.Tm(end)-S.Tm(1));
% implicit (backward) marching -> rectangle rule with end-of-step values
E_w   = -sum(dxv.*2*pi.*S.R.*S.qw_cons);
E_d   = sum(dxv.*S.diss);
S.energy_err = (E_b - E_w - E_d)/max(abs(E_w)+abs(E_d), eps);
S.Qheat_total = -E_w*p.Ntube;                            % [W], + = removed
S.Qt = Qt;
S.Re = p.rho*(Qt/(pi*p.R0^2))*2*p.R0/S.mufield(1,1);    % based on centreline mu at inlet
S.Re_wall = p.rho*(Qt/(pi*p.R0^2))*2*p.R0/S.mufield(end-1,1);
end

% ======================================================================
function [type, U, Text] = wall_bc(p, R, Tc)
% Returns BC type, the coefficient (U [W/m2/K] or q [W/m2]) and Text,
% all referred to the blood-side surface at radius R.
Rclot = 0;                                           % m^2 K / W
if R < p.R0*(1-1e-12)
    Rclot = R*log(p.R0/R)/p.clot.k;
end
switch lower(p.bc)
    case 't'
        if Rclot == 0, type = 'dirichlet'; U = Inf; else, type = 'robin'; U = 1/Rclot; end
        Text = p.T_wall;
    case 'q'
        type = 'neumann'; U = p.q_wall*p.R0/R; Text = NaN;   % same heat per length
    case 'coupled'
        Ro    = p.R0 + p.t_wall;
        Rtot  = Rclot + R*log(Ro/p.R0)/p.k_wall + R/(Ro*p.h_c);
        type  = 'robin_coupled'; U = 1/Rtot; Text = Tc;
    otherwise
        error('cp_solve:bc','Unknown bc "%s"',p.bc);
end
end

function [I,J,V,nz] = add(I,J,V,nz,ii,jj,vv)
n = numel(ii);
I(nz+1:nz+n) = ii; J(nz+1:nz+n) = jj; V(nz+1:nz+n) = vv; nz = nz+n;
end
