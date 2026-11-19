function V = cp_validate_transient()
%CP_VALIDATE_TRANSIENT  Verification of CP_SOLVE_TRANSIENT.
%   TV1  steady limit  == CP_SOLVE (all wall conditions); counter-current
%        convergence under grid refinement
%   TV2  exact discrete energy conservation (blood + wall + coolant)
%   TV3  pure-advection residence-time distribution of Poiseuille flow:
%        F(t) = 1 - (tau/2t)^2 ,  tau = L/U_mean  (analytical)
%   TV4  temporal order of accuracy (BDF2)
outdir = fullfile(fileparts(mfilename('fullpath')),'..','results');
if ~exist(outdir,'dir'), mkdir(outdir); end
fprintf('\n=========== TRANSIENT VALIDATION ===========\n');
p0 = cp_params(); p0.rheology = 'carreau-yasuda';
grd = struct('Nr',61,'Nx',200,'x_cluster',1.5);

%% TV1 steady limit
bcs = {'T','q','coupled-co','coupled-counter'}; V.tv1 = zeros(1,4);
for k = 1:4
    p = setbc(p0,bcs{k}); p.Nr = grd.Nr; p.Nx = grd.Nx; p.x_cluster = grd.x_cluster;
    tr = grd; tr.t_end = 0.5; tr.dt = 0.05;
    S  = cp_solve_transient(p,tr);  Ss = cp_solve(p);
    V.tv1(k) = max(abs(S.Tm_xt(:,1)'-Ss.Tm));
    fprintf('TV1 %-16s steady Tm: max |transient-code - cp_solve| = %.2e K, drift = %.1e K\n', ...
        bcs{k}, V.tv1(k), max(abs(S.Tm_xt(:,end)-S.Tm_xt(:,1))));
end
p = setbc(p0,'coupled-counter'); pf = p; pf.Nx = 4000; pf.Nr = 61; pf.x_cluster = 1.5;
Sf = cp_solve(pf);  Nxs = [100 200 400 800]; e = zeros(size(Nxs));
for k = 1:numel(Nxs)
    tr = struct('Nr',61,'Nx',Nxs(k),'x_cluster',1.5,'t_end',0.05,'dt',0.05);
    S = cp_solve_transient(p,tr);
    e(k) = abs(S.Tm_out(1) - Sf.Tm(end));
end
fprintf('    counter-current T_m(L) vs fine reference, Nx=%s : %s  (order %s)\n', ...
    mat2str(Nxs), mat2str(e,3), mat2str(log2(e(1:end-1)./e(2:end)),3));
V.tv1_counter = [Nxs; e];

%% TV2 energy conservation (backward Euler, coolant step)
for fa = {'counter','co'}
    p = setbc(p0,['coupled-' fa{1}]);
    tr = grd; tr.t_end = 20; tr.dt = 0.05; tr.scheme = 'BE'; tr.Tc_in = @(t) 28-8*(t>2);
    S = cp_solve_transient(p,tr);
    fprintf('TV2 %-8s energy residual %.2e J  (energy change %.2f J, rel %.1e)\n', fa{1}, ...
        S.energy_residual, S.E(end)-S.E(1), abs(S.energy_residual/(S.E(end)-S.E(1))));
end

%% TV3 residence-time distribution (pure advection, Newtonian)
p = p0; p.rheology = 'newtonian'; p.bc = 'q'; p.q_wall = 0; p.k = 1e-7; p.T_in = 37;
Um = p.Q_total/p.Ntube/(pi*p.R0^2);  tau = p.L/Um;
Nxs = [100 200 400 800]; err = zeros(size(Nxs)); fig = newfig([100 100 560 400]); hold on;
for k = 1:numel(Nxs)
    tr = struct('Nr',61,'Nx',Nxs(k),'x_cluster',1,'t_end',3*tau,'dt',tau/800,'T_in',@(t) 37+(t>0));
    S = cp_solve_transient(p,tr);
    F = S.Tm_out - 37;
    Fex = (S.t >= tau/2).*(1 - (tau./(2*max(S.t,eps))).^2);
    err(k) = trapz(S.t,abs(F-Fex))/tau;
    plot(S.t/tau,F,'DisplayName',sprintf('N_x = %d',Nxs(k)));
end
plot(S.t/tau,Fex,'k--','LineWidth',2,'DisplayName','exact 1-(\tau/2t)^2');
xlabel('t / \tau  (\tau = L/U_{mean})'); ylabel('F(t)'); grid on; legend('Location','southeast');
title('Residence-time distribution, Poiseuille flow');
exportgraphics(fig,fullfile(outdir,'VT3_residence_time.png'),'Resolution',200); close(fig);
fprintf('TV3 RTD L1 error, Nx=%s : %s (order %s)\n',mat2str(Nxs),mat2str(err,3), ...
    mat2str(log2(err(1:end-1)./err(2:end)),3));
V.tv3 = [Nxs; err];

%% TV4 temporal convergence (BDF2), smooth coolant ramp
p = setbc(p0,'coupled-counter');
ramp = @(t) 28 - 8*min(max(t-1,0),1).^2.*(3-2*min(max(t-1,0),1));   % smoothstep 1->2 s
dts = [0.2 0.1 0.05 0.025]; tq = 4;
tr = grd; tr.t_end = tq; tr.dt = 0.0025; tr.Tc_in = ramp;
Sref = cp_solve_transient(p,tr); yref = Sref.Tm_out(end);
et = zeros(size(dts));
for k = 1:numel(dts)
    tr.dt = dts(k); S = cp_solve_transient(p,tr); et(k) = abs(S.Tm_out(end)-yref);
end
fprintf('TV4 BDF2 error in T_m,out(t=%g s), dt=%s : %s (order %s)\n',tq,mat2str(dts), ...
    mat2str(et,3),mat2str(log2(et(1:end-1)./et(2:end)),3));
V.tv4 = [dts; et];
fig = newfig([100 100 480 380]); loglog(dts,et,'o-',dts,et(end)*(dts/dts(end)).^2,'k--'); grid on;
xlabel('\Delta t [s]'); ylabel('|error| [K]'); legend('BDF2','slope 2','Location','northwest');
title('Temporal convergence');
exportgraphics(fig,fullfile(outdir,'VT4_time_convergence.png'),'Resolution',200); close(fig);
fprintf('============================================\n');
end

function p = setbc(p,bc)
switch bc
    case 'T', p.bc = 'T';
    case 'q', p.bc = 'q';
    case 'coupled-counter', p.bc = 'coupled'; p.flow_arrangement = 'counter';
    case 'coupled-co',      p.bc = 'coupled'; p.flow_arrangement = 'co';
end
end

function f = newfig(pos)
f = figure('Visible','off','Position',pos,'Color','w'); try, theme(f,'light'); catch, end
end
