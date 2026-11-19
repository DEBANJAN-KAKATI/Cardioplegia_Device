function S = cp_solve_transient(p, tr)
%CP_SOLVE_TRANSIENT  Transient 2-D (r,x) conjugate heat transfer in the
%   cardioplegia heat exchanger: blood + tube wall + coolant, with
%   time-dependent operating inputs (digital-twin core).
%
%   Blood (per tube, mapped coordinate eta = r/R(x)):
%     rho cp ( T_t + u T_x - G_x/(R^2 eta) T_eta ) = (1/(R^2 eta)) (k eta T_eta)_eta + mu gd^2
%   Wall (lumped per unit length, only for bc = 'coupled'):
%     rho_w c_w A_w dTs/dt = 2 pi R [ U1 (T_wall - Ts) - U2 (Ts - Tc) ]
%   Coolant (1-D, co- or counter-current):
%     rho_c c_c A_c dTc/dt +/- mc c_c dTc/dx = 2 pi R U2 (Ts - Tc)
%
%   Momentum is quasi-steady (viscous time R^2/nu ~ 0.2 s << thermal times),
%   so the velocity field follows Q(t) instantaneously (lubrication model of
%   CP_VELOCITY).  Space: the same FV / upwind discretisation as CP_SOLVE
%   (so the steady limit reproduces CP_SOLVE).  Time: BDF2 (or backward
%   Euler), fully implicit, one sparse LU per change of the flow rates.
%   With mu(T) the viscosity is lagged by one time step (semi-implicit).
%
%   tr (all optional):
%     t_end, dt, scheme ('BDF2'|'BE'), Nr, Nx, x_cluster,
%     init ('steady' at t=0 inputs | 'uniform' with T0), T0, save_times,
%     T_in(t), Q_total(t), Tc_in(t), Qc_total(t), T_wall(t), q_wall(t)
%       -> function handles of time [s]; default: constants from p.
%
%   S: t, Tm_out, Tc_out, Tmax, Twall_max, p_in_mmHg, dp_mmHg, Qheat_W,
%      Tm_xt (Nx x Nt), snapshots, steady-init solution, energy residual.

if nargin < 2, tr = struct; end
tr = fill(tr, struct('t_end',30,'dt',0.05,'scheme','BDF2','Nr',61,'Nx',200, ...
    'x_cluster',1.5,'init','steady','T0',[],'save_times',[], ...
    'T_in',@(t) p.T_in+0*t,'Q_total',@(t) p.Q_total+0*t,'Tc_in',@(t) p.Tc_in+0*t, ...
    'Qc_total',@(t) p.Qc_total+0*t,'T_wall',@(t) p.T_wall+0*t,'q_wall',@(t) p.q_wall+0*t));

% ---------- grids
s   = linspace(0,1,tr.Nr)';
g.eta  = sin(pi/2*s);  N = numel(g.eta);
g.x    = p.L*linspace(0,1,tr.Nx+1).^tr.x_cluster;  Nx = numel(g.x);
g.h    = diff(g.eta);  g.etaf = 0.5*(g.eta(1:end-1)+g.eta(2:end));
g.V    = ([g.etaf;1].^2 - [0;g.etaf].^2)/2;
g.c    = [0 diff(g.x)];                       % cell length owned by column i
g.R    = cp_radius(g.x,p);
hm = g.h(1:end-1); hp = g.h(2:end);
g.Dm = -hp./(hm.*(hm+hp)); g.D0 = (hp-hm)./(hm.*hp); g.Dp = hm./(hp.*(hm+hp));
h1 = g.h(end); h2 = g.h(end-1);
g.Dw = [h1/(h2*(h1+h2)), -(h1+h2)/(h1*h2), (2*h1+h2)/(h1*(h1+h2))];
g.N = N; g.Nx = Nx; g.nb = N*Nx;
g.iTc = g.nb + (1:Nx); g.iTs = g.nb + Nx + (1:Nx); g.Ntot = g.nb + 2*Nx;
coupled = strcmpi(p.bc,'coupled');
counter = coupled && strcmpi(p.flow_arrangement,'counter');

% ---------- time grid
nt = round(tr.t_end/tr.dt);  t = (0:nt)*tr.dt;

% ---------- initial state
inp = inputs(tr,0);
Vel = columns(p,g,inp.Q_total/p.Ntube,[]);
[A,b,M] = assemble(p,g,Vel,inp,coupled,counter);
if strcmpi(tr.init,'steady')
    X = A\b;
    if p.Tdep
        for it = 1:50
            Vel = columns(p,g,inp.Q_total/p.Ntube,reshape(X(1:g.nb),N,Nx));
            [A,b,M] = assemble(p,g,Vel,inp,coupled,counter);
            Xn = A\b;  d = max(abs(Xn-X));  X = Xn;
            if d < 1e-8, break; end
        end
    end
else
    X = tr.T0*ones(g.Ntot,1);
    X(1:N) = inp.T_in;
end
S.X_init = X;

% ---------- storage
S.t = t; S.x = g.x; S.eta = g.eta; S.R = g.R;
z = zeros(1,nt+1);
S.Tm_out = z; S.Tc_out = z; S.Tmax = z; S.Twall_max = z; S.p_in_mmHg = z;
S.dp_mmHg = z; S.Qheat_W = z; S.tau_w_max = z; S.E = z; S.Ein = z;
S.Tm_xt = zeros(Nx,nt+1); S.Tc_xt = S.Tm_xt; S.inputs = repmat(inp,1,nt+1);
S.snap = struct('t',{},'T',{},'Tm',{},'Tc',{},'Ts',{});
record(1);

Xm1 = X;  Qprev = inp.Q_total; Qcprev = inp.Qc_total; coefPrev = NaN; dA = [];
for n = 1:nt
    inp = inputs(tr,t(n+1));
    newV = p.Tdep || abs(inp.Q_total-Qprev) > 1e-12*Qprev;
    if newV
        Tl = reshape(X(1:g.nb),N,Nx);             % lagged temperature
        Vel = columns(p,g,inp.Q_total/p.Ntube,Tl);
    end
    newA = newV || abs(inp.Qc_total-Qcprev) > 1e-12*max(Qcprev,eps);
    [A,b,M] = assemble(p,g,Vel,inp,coupled,counter);
    if strcmpi(tr.scheme,'BE') || n == 1
        coef = 1;   rhs = b + M.*X/tr.dt;
    else
        coef = 1.5; rhs = b + M.*(2*X - 0.5*Xm1)/tr.dt;
    end
    if newA || coef ~= coefPrev || isempty(dA)
        dA = decomposition(A + spdiags(coef*M/tr.dt,0,g.Ntot,g.Ntot),'lu');
        coefPrev = coef;
    end
    Xm1 = X;  X = dA\rhs;
    Qprev = inp.Q_total; Qcprev = inp.Qc_total;
    S.inputs(n+1) = inp;
    % discrete energy supplied during the step (for the conservation check)
    S.Ein(n+1) = S.Ein(n) + tr.dt*boundary_power(X);
    record(n+1);
end
S.energy_residual = (S.E(end)-S.E(1)) - S.Ein(end);   % exact (0) for BE, clean tube
S.p = p; S.tr = tr; S.g = g;

% ===================== nested helpers ===================================
    function record(k)
        T  = reshape(X(1:g.nb),N,Nx);
        Tc = X(g.iTc)'; Ts = X(g.iTs)';
        Tm = trapz(g.eta, Vel.u.*T.*g.eta) ./ trapz(g.eta, Vel.u.*g.eta);
        Tm(1) = T(1,1);
        S.Tm_xt(:,k) = Tm';  S.Tc_xt(:,k) = Tc';
        S.Tm_out(k) = Tm(end);
        if coupled, S.Tc_out(k) = ternary(counter, Tc(1), Tc(end)); else, S.Tc_out(k) = NaN; end
        S.Tmax(k) = max(T(:));  S.Twall_max(k) = max(T(N,2:end));
        dp = trapz(g.x, 2*Vel.tau_w./g.R);
        S.dp_mmHg(k) = dp/133.322;  S.p_in_mmHg(k) = p.p_out_mmHg + dp/133.322;
        S.tau_w_max(k) = max(Vel.tau_w);
        % heat removed from blood through the wall, all tubes [W]
        qw = wallflux(T,Ts);
        S.Qheat_W(k) = p.Ntube*sum(g.c.*2*pi.*g.R.*qw);
        % stored energy (per tube, discrete, relative to 0 C)
        Eb = sum(g.c.*2*pi.*g.R.^2.*(p.rho*p.cp*(g.V'*T)));
        Ew = 0; Ec = 0;
        if coupled
            Aw = pi*((p.R0+p.t_wall)^2 - p.R0^2);
            Ew = sum(g.c.*p.rho_w*p.c_w*Aw*p.wallCap.*Ts);
            if counter, cc = [g.c(2:end) 0]; else, cc = g.c; end
            Ec = sum(cc.*p.rho_c*p.cp_c*p.Ac.*Tc);
        end
        S.E(k) = Eb + Ew + Ec;
        if any(abs(t(k)-tr.save_times) < tr.dt/2)
            S.snap(end+1) = struct('t',t(k),'T',T,'Tm',Tm,'Tc',Tc,'Ts',Ts);
        end
    end

    function qw = wallflux(T,Ts)
        qw = zeros(1,Nx);
        for i = 2:Nx
            [type,U1,~,Text] = bcinfo(p,g.R(i),inp);
            switch type
                case 'dirichlet', qw(i) = -(p.k/g.R(i))*(g.Dw*T(N-2:N,i));
                case 'neumann',   qw(i) = U1;
                case 'robin',     qw(i) = U1*(T(N,i)-Text);
                case 'coupled',   qw(i) = U1*(T(N,i)-Ts(i));
            end
        end
    end

    function P = boundary_power(X)
        % net power entering the (blood+wall+coolant) system per tube [W]
        T  = reshape(X(1:g.nb),N,Nx);
        Ri = g.R(end);
        P  = 2*pi*Ri^2*p.rho*p.cp*sum(g.V.*Vel.u(:,end).*(T(:,1)-T(:,end)));
        if p.viscDiss
            P = P + sum(g.c.*2*pi.*g.R.^2.*(g.V'*Vel.Phi));
        end
        if coupled
            mcc = p.rho_c*p.cp_c*inp.Qc_total/p.Ntube;
            Tc = X(g.iTc);
            if counter, P = P + mcc*(Tc(end)-Tc(1)); else, P = P + mcc*(Tc(1)-Tc(end)); end
        else
            P = NaN;   % open wall: not a closed balance
        end
    end
end

% =======================================================================
function inp = inputs(tr,t)
inp = struct('T_in',tr.T_in(t),'Q_total',tr.Q_total(t),'Tc_in',tr.Tc_in(t), ...
    'Qc_total',tr.Qc_total(t),'T_wall',tr.T_wall(t),'q_wall',tr.q_wall(t));
end

function V = columns(p,g,Qt,T)
% quasi-steady velocity in every column
N = g.N; Nx = g.Nx;
V.u = zeros(N,Nx); V.Phi = V.u; V.tau_w = zeros(1,Nx);
if ~p.Tdep || isempty(T)
    [Ru,~,ic] = unique(g.R);
    for k = 1:numel(Ru)
        if p.Tdep, Tk = p.T_ref*ones(N,1); else, Tk = zeros(N,1); end
        Vk = cp_velocity(Ru(k),Qt,g.eta,Tk,p);
        cols = (ic==k);
        V.u(:,cols) = repmat(Vk.u,1,nnz(cols));
        V.Phi(:,cols) = repmat(Vk.Phi,1,nnz(cols));
        V.tau_w(cols) = Vk.tau_w;
    end
else
    for i = 1:Nx
        Vk = cp_velocity(g.R(i),Qt,g.eta,T(:,i),p);
        V.u(:,i) = Vk.u; V.Phi(:,i) = Vk.Phi; V.tau_w(i) = Vk.tau_w;
    end
end
G = g.R.^2 .* cumtrapz(g.eta, V.u.*g.eta);
V.W = zeros(N,Nx);
for i = 2:Nx
    V.W(2:N-1,i) = (G(2:N-1,i)-G(2:N-1,i-1))/g.c(i) ./ (g.R(i)^2*g.eta(2:N-1));
end
end

function [type,U1,U2,Text] = bcinfo(p,R,inp)
Rclot = 0; if R < p.R0*(1-1e-12), Rclot = R*log(p.R0/R)/p.clot.k; end
U2 = NaN; Text = NaN;
switch lower(p.bc)
    case 't'
        Text = inp.T_wall;
        if Rclot == 0, type = 'dirichlet'; U1 = Inf; else, type = 'robin'; U1 = 1/Rclot; end
    case 'q'
        type = 'neumann'; U1 = inp.q_wall*p.R0/R;
    case 'coupled'
        Ro = p.R0 + p.t_wall; Rm = 0.5*(p.R0+Ro);
        U1 = 1/(Rclot + R*log(Rm/p.R0)/p.k_wall);
        U2 = 1/(R*log(Ro/Rm)/p.k_wall + R/(Ro*p.h_c));
        type = 'coupled';
end
end

function [A,b,M] = assemble(p,g,V,inp,coupled,counter)
N = g.N; Nx = g.Nx; Ntot = g.Ntot; rc = p.rho*p.cp;
nzmax = Nx*(N*12+20); I = zeros(nzmax,1); J = I; Vv = I; nz = 0;
b = zeros(Ntot,1); M = zeros(Ntot,1); isDir = false(Ntot,1);
    function ad(ii,jj,vv)
        n_ = numel(ii); I(nz+1:nz+n_) = ii; J(nz+1:nz+n_) = jj; Vv(nz+1:nz+n_) = vv; nz = nz+n_;
    end
% inlet column (Dirichlet)
isDir(1:N) = true; b(1:N) = inp.T_in;
jj = (2:N-1)'; f = (1:N-1)';
for i = 2:Nx
    R = g.R(i); dx = g.c(i); u = V.u(:,i); W = V.W(:,i);
    k = (i-1)*N + (1:N)';
    cd = rc*g.V.*u/dx;
    ad(k,k,cd); ad(k,k-N,-cd);
    M(k) = rc*g.V;
    cw = -rc*g.V(jj).*W(jj);
    ad(k(jj),k(jj)-1,cw.*g.Dm); ad(k(jj),k(jj),cw.*g.D0); ad(k(jj),k(jj)+1,cw.*g.Dp);
    a = p.k*g.etaf./g.h/R^2;
    ad(k(f),k(f),a); ad(k(f),k(f)+1,-a); ad(k(f)+1,k(f)+1,a); ad(k(f)+1,k(f),-a);
    if p.viscDiss, b(k) = b(k) + g.V.*V.Phi(:,i); end
    kN = k(N);
    [type,U1,U2,Text] = bcinfo(p,R,inp);
    switch type
        case 'dirichlet', isDir(kN) = true; b(kN) = Text; M(kN) = 0;
        case 'neumann',   b(kN) = b(kN) - U1/R;
        case 'robin',     ad(kN,kN,U1/R); b(kN) = b(kN) + U1/R*Text;
        case 'coupled'
            % coolant cell [x(i-1),x(i)] is stored at its OUTFLOW end:
            % node i (co-current) or node i-1 (counter-current)  -> exact
            % discrete energy conservation with first-order upwinding
            its = g.iTs(i); c = g.c(i);
            if counter, itc = g.iTc(i-1); else, itc = g.iTc(i); end
            ad(kN,kN,U1/R); ad(kN,its,-U1/R);
            % wall node (row weighted by cell length c)
            Aw = pi*((p.R0+p.t_wall)^2 - p.R0^2);
            M(its) = c*p.rho_w*p.c_w*Aw*p.wallCap;
            ad(its,its,c*2*pi*R*(U1+U2)); ad(its,kN,-c*2*pi*R*U1); ad(its,itc,-c*2*pi*R*U2);
            % coolant node
            M(itc) = c*p.rho_c*p.cp_c*p.Ac;
            ad(itc,itc,c*2*pi*R*U2); ad(itc,its,-c*2*pi*R*U2);
    end
end
% coolant convection / inlets, wall node at x=0
if coupled
    mcc = p.rho_c*p.cp_c*inp.Qc_total/p.Ntube;
    if counter
        for i = 1:Nx-1
            ad(g.iTc(i),g.iTc(i),mcc); ad(g.iTc(i),g.iTc(i+1),-mcc);
        end
        isDir(g.iTc(Nx)) = true; b(g.iTc(Nx)) = inp.Tc_in; M(g.iTc(Nx)) = 0;
    else
        for i = 2:Nx
            ad(g.iTc(i),g.iTc(i),mcc); ad(g.iTc(i),g.iTc(i-1),-mcc);
        end
        isDir(g.iTc(1)) = true; b(g.iTc(1)) = inp.Tc_in;
    end
    ad(g.iTs(1),g.iTs(1),1); ad(g.iTs(1),g.iTs(2),-1);          % Ts(0) = Ts(dx)
else
    isDir([g.iTc g.iTs]) = true;                               % unused unknowns = 0
end
A = sparse(I(1:nz),J(1:nz),Vv(1:nz),Ntot,Ntot);
Z = spdiags(double(~isDir),0,Ntot,Ntot);
A = Z*A + spdiags(double(isDir),0,Ntot,Ntot);
M(isDir) = 0;
end

function s = fill(s,d)
f = fieldnames(d); for i = 1:numel(f), if ~isfield(s,f{i}), s.(f{i}) = d.(f{i}); end, end
end
function o = ternary(c,a,b), if c, o = a; else, o = b; end, end
