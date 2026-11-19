function V = cp_validate()
%CP_VALIDATE  Verification of CP_SOLVE against analytical solutions.
%   V1  velocity profiles : Hagen-Poiseuille, power-law, CY->Newtonian limit
%   V2  Graetz problem (const. wall T)   : series solution, Nu_inf = 3.657
%   V3  const. heat flux                  : Tm linear, Nu_inf = 48/11
%   V4  power-law const. flux             : Nu_inf = 8(5n+1)(3n+1)/(31n^2+12n+1)
%   V5  adiabatic wall + dissipation      : Tm(L)-Tin = dp/(rho cp)
%   V6  grid convergence
outdir = fullfile(fileparts(mfilename('fullpath')),'..','results');
if ~exist(outdir,'dir'), mkdir(outdir); end
p0 = cp_params();
R  = p0.R0;  Qt = p0.Q_total/p0.Ntube;  Um = Qt/(pi*R^2);
fprintf('\n================ VALIDATION ================\n');

%% V1 velocity
eta = sin(pi/2*linspace(0,1,201))';
p = p0; p.rheology = 'newtonian';
Vn = cp_velocity(R,Qt,eta,37*ones(size(eta)),p);
e1 = max(abs(Vn.u - 2*Um*(1-eta.^2)))/(2*Um);
e1t = abs(Vn.tau_w - 4*p.mu_N*Qt/(pi*R^3))/(4*p.mu_N*Qt/(pi*R^3));
p.rheology = 'powerlaw'; n = p.PL.n;
Vp = cp_velocity(R,Qt,eta,37*ones(size(eta)),p);
up = Um*(3*n+1)/(n+1)*(1-eta.^((n+1)/n));
e2 = max(abs(Vp.u-up))/max(up);
twp = p.PL.K*((3*n+1)/n*Qt/(pi*R^3))^n;
e2t = abs(Vp.tau_w-twp)/twp;
p.rheology = 'carreau-yasuda'; p.CY.mu0 = p.mu_N; p.CY.muinf = p.mu_N;
Vc = cp_velocity(R,Qt,eta,37*ones(size(eta)),p);
e3 = max(abs(Vc.u - 2*Um*(1-eta.^2)))/(2*Um);
fprintf('V1  Newtonian u err %.2e, tau_w err %.2e\n',e1,e1t);
fprintf('    Power-law u err %.2e, tau_w err %.2e\n',e2,e2t);
fprintf('    CY (mu0=muinf) vs Poiseuille u err %.2e\n',e3);
V.v1 = [e1 e1t e2 e2t e3];

%% V2 Graetz (constant wall temperature), Newtonian
p = p0; p.rheology='newtonian'; p.bc='T'; p.T_in=42; p.T_wall=32; p.L=0.3;
p.Nx = 800; p.Nr = 121;
S = cp_solve(p);
alpha = p.k/(p.rho*p.cp);
zeta  = alpha*S.x/(2*Um*R^2);
[thm, Nug] = graetz_series(zeta, 80);
thS = (S.Tm - p.T_wall)/(p.T_in - p.T_wall);
xp  = S.x/(2*R*Um*2*R/alpha);                       % x+ = x/(D Pe)
mask = zeta > 1e-3;
e4 = max(abs(thS(mask)-thm(mask)));
e5 = max(abs(S.Nu(mask)-Nug(mask))./Nug(mask));
fprintf('V2  Graetz: max |theta_m err| = %.2e, max rel Nu err (zeta>1e-3) = %.2e\n',e4,e5);
fprintf('    Nu(end) = %.4f   (exact 3.6568)\n',S.Nu(end));
V.v2 = [e4 e5 S.Nu(end)];
fig = figure('Visible','off','Position',[100 100 1000 380]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile; semilogx(xp(2:end),thm(2:end),'k-','LineWidth',1.5); hold on;
semilogx(xp(2:8:end),thS(2:8:end),'ro','MarkerSize',4);
xlabel('x^+ = x/(D Pe)'); ylabel('\theta_m'); legend('Graetz series','present FV','Location','southwest'); grid on;
title('Cup-mixing temperature');
nexttile; loglog(xp(2:end),Nug(2:end),'k-','LineWidth',1.5); hold on;
loglog(xp(2:8:end),S.Nu(2:8:end),'ro','MarkerSize',4); yline(3.6568,'b--');
xlabel('x^+'); ylabel('Nu_x'); legend('Graetz series','present FV','Nu_\infty = 3.657'); grid on;
title('Local Nusselt number (constant T_w)');
exportgraphics(fig, fullfile(outdir,'V2_graetz_validation.png'),'Resolution',200); close(fig);

%% V3 constant flux, Newtonian
p = p0; p.rheology='newtonian'; p.bc='q'; p.L=0.3; p.Nx=800; p.Nr=121;
S = cp_solve(p);
Tm_ex = p.T_in - p.q_wall*2*pi*R*S.x/(p.rho*p.cp*Qt);
e6 = max(abs(S.Tm - Tm_ex));
fprintf('V3  Const flux: max |Tm - exact| = %.2e K, Nu(end) = %.4f (exact %.4f)\n',e6,S.Nu(end),48/11);
V.v3 = [e6 S.Nu(end)];

%% V4 power-law, constant flux: Nu_inf analytical
p.rheology='powerlaw';
nl = [0.5 0.708 1.0]; Nnum = zeros(size(nl)); Nex = Nnum;
for k = 1:numel(nl)
    p.PL.n = nl(k); p.PL.K = 0.017*(100)^(1-nl(k))*(p.mu_N/0.017);   % similar mu at 100 1/s
    S = cp_solve(p);
    Nnum(k) = S.Nu(end);
    n = nl(k); Nex(k) = 8*(5*n+1)*(3*n+1)/(31*n^2+12*n+1);
    fprintf('V4  Power-law n=%.3f : Nu_inf = %.4f  (analytical %.4f, err %.2e)\n',n,Nnum(k),Nex(k),abs(Nnum(k)-Nex(k))/Nex(k));
end
V.v4 = [Nnum; Nex];

%% V5 adiabatic wall with viscous dissipation:  rho cp Q dTm = Q dp
p = p0; p.rheology='carreau-yasuda'; p.bc='q'; p.q_wall=0; p.viscDiss=true; p.T_in=37;
S = cp_solve(p);
dT_num = S.Tm(end)-p.T_in;  dT_ex = S.dp/(p.rho*p.cp);
fprintf('V5  Adiabatic + dissipation (CY): dTm = %.4e K, dp/(rho cp) = %.4e K, rel err %.2e\n',dT_num,dT_ex,abs(dT_num-dT_ex)/dT_ex);
V.v5 = [dT_num dT_ex];

%% V6 grid convergence (Graetz, Nu at x+ = 0.01 and Tm(L))
p = p0; p.rheology='newtonian'; p.bc='T'; p.T_in=42; p.T_wall=32;
xq = 0.01*(2*R)*(2*R*Um/alpha);
Nrs = [31 61 121 241]; Nxs = [125 250 500 1000];
ErrNu = zeros(size(Nrs)); ErrT = ErrNu;
[thm_q, Nu_q] = graetz_series(alpha*xq/(2*Um*R^2), 80);
for k = 1:numel(Nrs)
    p.Nr = Nrs(k); p.Nx = Nxs(k);
    S = cp_solve(p);
    Nk = interp1(S.x(2:end),S.Nu(2:end),xq);
    Tk = (interp1(S.x,S.Tm,xq)-p.T_wall)/(p.T_in-p.T_wall);
    ErrNu(k) = abs(Nk-Nu_q)/Nu_q;  ErrT(k) = abs(Tk-thm_q);
    fprintf('V6  Nr=%4d Nx=%5d : rel err Nu = %.3e, err theta_m = %.3e, energy err = %.2e\n',Nrs(k),Nxs(k),ErrNu(k),ErrT(k),S.energy_err);
end
ordNu = log(ErrNu(1:end-1)./ErrNu(2:end))/log(2);
ordT  = log(ErrT(1:end-1)./ErrT(2:end))/log(2);
fprintf('    observed order (Nu): %s\n    observed order (Tm): %s\n',mat2str(ordNu,3),mat2str(ordT,3));
V.v6 = [Nrs; Nxs; ErrNu; ErrT];
fig = figure('Visible','off','Position',[100 100 500 380]);
loglog(Nrs,ErrNu,'o-',Nrs,ErrT,'s-','LineWidth',1.4); grid on;
xlabel('N_r  (N_x = 4.1 N_r)'); ylabel('error'); legend('rel. error Nu','abs. error \theta_m');
title('Grid convergence (Graetz, x^+=0.01)');
exportgraphics(fig, fullfile(outdir,'V6_grid_convergence.png'),'Resolution',200); close(fig);
fprintf('============================================\n');
end

% ----------------------------------------------------------------------
function [thm, Nu] = graetz_series(zeta, nm)
% Independent eigen-solution of  (eta phi')' + lam^2 eta (1-eta^2) phi = 0,
% phi'(0)=0, phi(1)=0  on a fine uniform cell-centred grid.
M = 4000; h = 1/M; e = ((1:M)'-0.5)*h; ef = (0:M)'*h;
main = (ef(1:end-1)+ef(2:end))/h^2;  main(end) = (ef(end-1) + 2*ef(end))/h^2;
off  = -ef(2:end-1)/h^2;
K = spdiags([[off;0] main [0;off]],[-1 0 1],M,M);
B = spdiags(e.*(1-e.^2),0,M,M);
[Phi,L] = eigs(K,B,nm,'smallestabs');
lam2 = diag(L); [lam2,o] = sort(lam2); Phi = Phi(:,o);
w  = e.*(1-e.^2)*h;
Sn = (w'*Phi)';                                     % int eta(1-eta^2) phi
Cn = Sn ./ ((w'*(Phi.^2))');                        % expansion of theta=1
E  = exp(-lam2*zeta(:)');                           % nm x Nz
num = (Cn.*Sn)'*E;  thm = 4*num;
Nu  = ((Cn.*Sn.*lam2)'*E)./(2*num);
thm = thm(:)'; Nu = Nu(:)';
end
