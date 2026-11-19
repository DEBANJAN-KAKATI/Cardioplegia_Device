function R = run_transient_cases()
%RUN_TRANSIENT_CASES  Digital-twin scenarios with the transient model.
%   D1  coolant inlet temperature step  (28 -> 20 C)         + FOPDT model
%   D2  blood inlet temperature step    (42 -> 43 C)         + FOPDT model
%   D3  blood flow step                 (300 -> 450 mL/min)
%   D4  REWARMING pump stop: blood flow 300 -> 15 mL/min for 30 s, then
%       restart -> time to exceed 40 C and the hot bolus sent to the patient
%   D5  pulsatile roller-pump flow      Q = Q0 (1 + 0.4 sin 2 pi f t), f = 1 Hz
%   Base configuration: cooling scenario, Carreau-Yasuda blood,
%   counter-current coupled water circuit, wall heat capacity included.
here = fileparts(mfilename('fullpath'));
outdir = fullfile(here,'..','results');
if ~exist(outdir,'dir'), mkdir(outdir); end
diary(fullfile(outdir,'transient_log.txt')); diary on; t0 = tic;
mm = 1e3; C = lines(6);

p = cp_params(); p.rheology = 'carreau-yasuda';
p.bc = 'coupled'; p.flow_arrangement = 'counter';
grd = struct('Nr',61,'Nx',200,'x_cluster',1.5,'dt',0.05);
Q0 = p.Q_total;  Um = Q0/p.Ntube/(pi*p.R0^2);
fprintf('Mean blood residence time L/U = %.2f s, radial diffusion time R^2/alpha = %.2f s\n', ...
    p.L/Um, p.R0^2/(p.k/(p.rho*p.cp)));

%% D1 coolant temperature step
tr = grd; tr.t_end = 40; tr.Tc_in = @(t) 28 - 8*(t>=2); tr.save_times = [2 3 4 6 10 40];
D1 = cp_solve_transient(p,tr);
F1 = fopdt(D1.t, D1.Tm_out, 2, -8);
fprintf('D1 Tc_in step -8 K : Tm_out %.2f -> %.2f C, FOPDT K=%.3f, dead time %.2f s, tau %.2f s, t95 = %.1f s\n', ...
    D1.Tm_out(1), D1.Tm_out(end), F1.K, F1.theta, F1.tau, F1.t95);

%% D2 blood inlet temperature step
tr = grd; tr.t_end = 40; tr.T_in = @(t) 42 + (t>=2);
D2 = cp_solve_transient(p,tr);
F2 = fopdt(D2.t, D2.Tm_out, 2, 1);
fprintf('D2 T_in step +1 K  : FOPDT K=%.3f, dead time %.2f s, tau %.2f s, t95 = %.1f s\n',F2.K,F2.theta,F2.tau,F2.t95);

%% D3 blood flow step
tr = grd; tr.t_end = 40; tr.Q_total = @(t) Q0*(1 + 0.5*(t>=2));
D3 = cp_solve_transient(p,tr);
F3 = fopdt(D3.t, D3.Tm_out, 2, 150);
fprintf('D3 flow step +50%%  : Tm_out %.2f -> %.2f C, p_in %.1f -> %.1f mmHg, tau %.2f s\n', ...
    D3.Tm_out(1),D3.Tm_out(end),D3.p_in_mmHg(1),D3.p_in_mmHg(end),F3.tau);

%% D4 rewarming with pump stop / restart
pr = p; pr.T_in = 32; pr.Tc_in = 42;
tr = grd; tr.t_end = 70; tr.Q_total = @(t) Q0*(1 - 0.95*(t>=5 & t<35));
tr.save_times = [5 15 35 36 40];
D4 = cp_solve_transient(pr,tr);
x40 = cp_crossing(D4.t, D4.Tm_out, pr.T_safe);
above = D4.Tm_out > pr.T_safe;
fprintf(['D4 pump stop (rewarming): Tm_out %.2f C before; reaches 40 C at t = %.1f s ' ...
    '(%.1f s after stop); time above 40 C = %.1f s; peak Tm_out = %.2f C; peak local T = %.2f C\n'], ...
    D4.Tm_out(1), x40, x40-5, sum(above)*grd.dt, max(D4.Tm_out), max(D4.Tmax));
% bolus after restart: blood volume delivered above 40 C
Qt_hist = arrayfun(@(s) s.Q_total, D4.inputs);
Vhot = trapz(D4.t, Qt_hist.*above)*1e6;
fprintf('    volume delivered above 40 C: %.1f mL\n', Vhot);

%% D5 pulsatile roller pump
tr = grd; tr.t_end = 12; tr.dt = 0.02; f = 1.0;
tr.Q_total = @(t) Q0*(1 + 0.4*sin(2*pi*f*t));
D5 = cp_solve_transient(p,tr);
last = D5.t > tr.t_end - 3;
Wo = p.R0*sqrt(2*pi*f*p.rho/p.CY.muinf);
fprintf('D5 pulsatile (f=%g Hz, Womersley %.2f): Tm_out %.3f +/- %.3f C, p_in %.1f +/- %.1f mmHg\n', ...
    f, Wo, mean(D5.Tm_out(last)), (max(D5.Tm_out(last))-min(D5.Tm_out(last)))/2, ...
    mean(D5.p_in_mmHg(last)), (max(D5.p_in_mmHg(last))-min(D5.p_in_mmHg(last)))/2);

%% ---------------- figures --------------------------------------------
% FT1: step responses + FOPDT fits
fg = newfig([50 50 1250 380]); tiledlayout(1,3,'TileSpacing','compact');
nexttile; plot(D1.t,D1.Tm_out,'Color',C(1,:)); hold on; plot(D1.t,F1.model,'k--');
plot(D1.t,D1.Tc_out,'Color',C(4,:)); grid on; xlabel('t [s]'); ylabel('T [^oC]');
legend('T_{m,out} (2-D model)','FOPDT fit','T_{c,out}'); title(sprintf('D1: T_{c,in} 28\\rightarrow20 ^oC  (\\tau=%.1f s, \\theta=%.1f s)',F1.tau,F1.theta));
nexttile; plot(D2.t,D2.Tm_out,'Color',C(2,:)); hold on; plot(D2.t,F2.model,'k--'); grid on;
xlabel('t [s]'); ylabel('T_{m,out} [^oC]'); legend('2-D model','FOPDT fit','Location','southeast');
title(sprintf('D2: T_{in} 42\\rightarrow43 ^oC  (K=%.2f, \\tau=%.1f s)',F2.K,F2.tau));
nexttile; X = D1.x*mm; contourf(D1.t, X, D1.Tm_xt, 30, 'LineColor','none'); colorbar;
xlabel('t [s]'); ylabel('x [mm]'); title('D1: T_m(x,t) [^oC]');
exportgraphics(fg,fullfile(outdir,'FT1_step_responses.png'),'Resolution',200); close(fg);

% FT2: flow step
fg = newfig([50 50 900 360]); tiledlayout(1,2,'TileSpacing','compact');
nexttile; plot(D3.t,D3.Tm_out,D3.t,D3.Tc_out); grid on; xlabel('t [s]'); ylabel('T [^oC]');
legend('T_{m,out}','T_{c,out}'); title('D3: blood flow 300\rightarrow450 mL/min');
nexttile; plot(D3.t,D3.p_in_mmHg); grid on; xlabel('t [s]'); ylabel('p_{in} [mmHg]'); title('Inlet pressure');
exportgraphics(fg,fullfile(outdir,'FT2_flow_step.png'),'Resolution',200); close(fg);

% FT3: pump stop during rewarming
fg = newfig([50 50 1250 650]); tiledlayout(2,3,'TileSpacing','compact');
nexttile([1 2]); yyaxis left; plot(D4.t,D4.Tm_out,D4.t,D4.Tmax,'-.'); ylabel('T [^oC]');
yline(40,'r-','40 ^oC'); yline(42,'m:','42 ^oC');
yyaxis right; plot(D4.t,Qt_hist*6e7,'k'); ylabel('Q_{blood} [mL/min]'); grid on; xlabel('t [s]');
legend('T_{m,out}','max local T','','','Q','Location','west');
title(sprintf('D4: rewarming, pump stop 5-35 s: T_{m,out}>40 ^oC for %.1f s, hot volume %.1f mL',sum(above)*grd.dt,Vhot));
nexttile; plot(D4.x*mm,vertcat(D4.snap.Tm)'); grid on; yline(40,'r:'); xlabel('x [mm]'); ylabel('T_m [^oC]');
legend(arrayfun(@(s) sprintf('t=%g s',s.t),D4.snap,'UniformOutput',false),'Location','southeast'); title('T_m(x) snapshots');
for k = [2 3 5]
    nexttile; s = D4.snap(k); contourf(D4.x*mm, D4.eta*p.R0*mm, s.T, 30,'LineColor','none'); colorbar;
    clim([32 42]); xlabel('x [mm]'); ylabel('r [mm]'); title(sprintf('T(r,x) at t = %g s',s.t));
end
exportgraphics(fg,fullfile(outdir,'FT3_pump_stop_rewarming.png'),'Resolution',200); close(fg);

% FT4: pulsatile
fg = newfig([50 50 1000 360]); tiledlayout(1,2,'TileSpacing','compact');
Qp = arrayfun(@(s) s.Q_total, D5.inputs)*6e7;
nexttile; yyaxis left; plot(D5.t,D5.Tm_out); ylabel('T_{m,out} [^oC]'); yyaxis right; plot(D5.t,Qp); ylabel('Q [mL/min]');
grid on; xlabel('t [s]'); title('D5: pulsatile roller-pump flow (1 Hz, \pm40 %)');
nexttile; plot(D5.t,D5.p_in_mmHg); grid on; xlabel('t [s]'); ylabel('p_{in} [mmHg]'); title('Inlet pressure');
exportgraphics(fg,fullfile(outdir,'FT4_pulsatile.png'),'Resolution',200); close(fg);

R = struct('D1',F1,'D2',F2,'D3',F3,'D4_t40',x40,'D4_time_above',sum(above)*grd.dt,'D4_Vhot_mL',Vhot);
save(fullfile(outdir,'transient_results.mat'),'R','D1','D2','D3','D4','D5','-v7.3');
fprintf('Transient cases done in %.1f s\n',toc(t0)); diary off;
end

% ----------------------------------------------------------------------
function F = fopdt(t,y,tstep,du)
% First-order-plus-dead-time fit  y = y0 + K du (1-exp(-(t-ts-theta)/tau))
y0 = y(find(t<tstep,1,'last'));  yinf = y(end);
K = (yinf-y0)/du;  m = t >= tstep;
model = @(q,tt) y0 + K*du*(1-exp(-max(tt-tstep-q(1),0)/q(2)));
cost  = @(q) sum((model(abs(q),t(m))-y(m)).^2);
q = abs(fminsearch(cost,[0.5 3],optimset('TolX',1e-6,'TolFun',1e-12,'MaxFunEvals',4000)));
F.K = K; F.theta = q(1); F.tau = q(2); F.model = model(q,t);
F.t95 = cp_crossing(t, abs(y-y0), 0.95*abs(yinf-y0)) - tstep;
end

function f = newfig(pos)
f = figure('Visible','off','Position',pos,'Color','w'); try, theme(f,'light'); catch, end
end
