function make_figures(RES, outdir, base, sweep, Qc_list, sev, dl, dis, Qlist, dmap, Qmap)
%MAKE_FIGURES  Publication figures for cases 4(a)-(i).
set(groot,'defaultAxesFontSize',10,'defaultLineLineWidth',1.4, ...
    'defaultAxesTickLabelInterpreter','tex','defaultTextInterpreter','tex');
mm = 1e3;  C = lines(6);
sv = @(f,name) exportgraphics(f, fullfile(outdir,name),'Resolution',200);

%% Fig 1  (a),(b) Newtonian: constant wall T vs constant flux
A = RES.cool_newtonian_T;  B = RES.cool_newtonian_q;
f = newfig([50 50 1200 360]); tiledlayout(1,3,'TileSpacing','compact');
nexttile; plot(A.x*mm,A.Tm,'Color',C(1,:)); hold on; plot(A.x*mm,A.Twall,'--','Color',C(1,:));
plot(B.x*mm,B.Tm,'Color',C(2,:)); plot(B.x*mm,B.Twall,'--','Color',C(2,:));
yline(40,'r:'); yline(37,'k:'); xlabel('x [mm]'); ylabel('T [^oC]');
legend('(a) T_m','(a) T_{wall}','(b) T_m','(b) T_{wall}','40 ^oC','37 ^oC','Location','northeast'); grid on;
title('Cup-mixing & wall temperature');
nexttile; semilogy(A.x(2:end)*mm,A.Nu(2:end)); hold on; semilogy(B.x(2:end)*mm,B.Nu(2:end));
yline(3.657,'k:'); yline(4.364,'k--'); ylim([3 30]); grid on;
xlabel('x [mm]'); ylabel('Nu_x'); legend('(a) const T_w','(b) const q_w','3.657','4.364'); title('Local Nusselt number');
nexttile; contourf(A.x*mm, A.eta*base.R0*mm, A.Tfield, 30,'LineColor','none'); colorbar;
xlabel('x [mm]'); ylabel('r [mm]'); title('(a) T(r,x) [^oC]');
sv(f,'F1_cases_a_b_newtonian.png'); close(f);

%% Fig 2  (c) rheology
N = RES.cool_newtonian_T; Cy = RES.cool_carreau_yasuda_T; PL = RES.cool_powerlaw_T;
f = newfig([50 50 1200 700]); tiledlayout(2,3,'TileSpacing','compact');
nexttile; Um = N.Qt/(pi*base.R0^2);
plot(N.ufield(:,end)/Um, N.eta,'Color',C(1,:)); hold on;
plot(Cy.ufield(:,end)/Um, Cy.eta,'Color',C(2,:)); plot(PL.ufield(:,end)/Um, PL.eta,'Color',C(3,:));
xlabel('u/U_{mean}'); ylabel('r/R'); legend('Newtonian','Carreau-Yasuda','power law'); grid on; title('Velocity profile');
nexttile; semilogx(N.mufield(:,end)*1e3,N.eta,'Color',C(1,:)); hold on;
semilogx(Cy.mufield(:,end)*1e3,Cy.eta,'Color',C(2,:)); semilogx(PL.mufield(:,end)*1e3,PL.eta,'Color',C(3,:));
xlabel('\mu [mPa s]'); ylabel('r/R'); grid on; title('Apparent viscosity');
nexttile; plot(N.x*mm,N.Tm,'Color',C(1,:)); hold on; plot(Cy.x*mm,Cy.Tm,'Color',C(2,:));
plot(RES.cool_newtonian_q.x*mm,RES.cool_newtonian_q.Tm,'--','Color',C(1,:));
plot(RES.cool_carreau_yasuda_q.x*mm,RES.cool_carreau_yasuda_q.Tm,'--','Color',C(2,:));
yline(40,'r:'); yline(37,'k:'); xlabel('x [mm]'); ylabel('T_m [^oC]'); grid on;
legend('N, const T_w','CY, const T_w','N, const q_w','CY, const q_w'); title('Cup-mixing temperature');
nexttile; semilogy(N.x(2:end)*mm,N.Nu(2:end),'Color',C(1,:)); hold on;
semilogy(Cy.x(2:end)*mm,Cy.Nu(2:end),'Color',C(2,:));
semilogy(RES.cool_newtonian_q.x(2:end)*mm,RES.cool_newtonian_q.Nu(2:end),'--','Color',C(1,:));
semilogy(RES.cool_carreau_yasuda_q.x(2:end)*mm,RES.cool_carreau_yasuda_q.Nu(2:end),'--','Color',C(2,:));
ylim([3 30]); grid on; xlabel('x [mm]'); ylabel('Nu_x'); title('Nusselt number');
nexttile; plot(N.x*mm, (Cy.Tm-N.Tm),'Color',C(2,:)); hold on;
plot(N.x*mm, RES.cool_carreau_yasuda_q.Tm-RES.cool_newtonian_q.Tm,'--','Color',C(2,:));
xlabel('x [mm]'); ylabel('T_m^{CY}-T_m^{N} [K]'); legend('const T_w','const q_w'); grid on; title('Non-Newtonian effect on T_m');
nexttile; plot(N.x*mm,N.tau_w,'Color',C(1,:)); hold on; plot(Cy.x*mm,Cy.tau_w,'Color',C(2,:)); plot(PL.x*mm,PL.tau_w,'Color',C(3,:));
xlabel('x [mm]'); ylabel('\tau_w [Pa]'); grid on; title('Wall shear stress');
sv(f,'F2_case_c_rheology.png'); close(f);

%% Fig 3  (d) coupled coolant circuit
f = newfig([50 50 1200 380]); tiledlayout(1,3,'TileSpacing','compact');
ks = {'cool_newtonian_coupled_counter','cool_carreau_yasuda_coupled_counter','cool_newtonian_coupled_co','cool_carreau_yasuda_coupled_co'};
lb = {'N counter','CY counter','N co','CY co'}; ls = {'-','-','--','--'}; cc = [1 2 1 2];
nexttile; hold on;
for k=1:4, S=RES.(ks{k}); plot(S.x*mm,S.Tm,ls{k},'Color',C(cc(k),:)); end
for k=[1 3], S=RES.(ks{k}); plot(S.x*mm,S.Tc,ls{k},'Color',C(4,:)); end
yline(40,'r:'); yline(37,'k:'); xlabel('x [mm]'); ylabel('T [^oC]'); grid on;
legend([lb {'T_c counter','T_c co'}],'Location','northeast'); title('Blood T_m and coolant T_c');
nexttile; hold on;
for k=1:4, S=RES.(ks{k}); plot(S.x*mm,S.Twall,ls{k},'Color',C(cc(k),:)); end
xlabel('x [mm]'); ylabel('T_{wall} [^oC]'); grid on; legend(lb); title('Blood-side wall temperature T_w(x)');
nexttile; semilogx(Qc_list,sweep(:,1),'o-',Qc_list,sweep(:,2),'s--'); grid on; yline(37,'k:');
xlabel('coolant flow [L/min]'); ylabel('T_m(L) [^oC]'); legend('counter-current','co-current'); title('Outlet temperature vs coolant flow (CY)');
sv(f,'F3_case_d_coupled.png'); close(f);

%% Fig 4  (e) cup-mixing temperature, length limits
f = newfig([50 50 1100 400]); tiledlayout(1,2,'TileSpacing','compact');
bcs = {'T','q','coupled_counter','coupled_co'}; rh = {'newtonian','carreau_yasuda'};
for sc = {'cool','rewarm'}
    nexttile; hold on; k=0; hh = gobjects(8,1);
    for b=1:4, for r=1:2
        k=k+1; S = RES.([sc{1} '_' rh{r} '_' bcs{b}]);
        hh(k) = plot(S.x*mm,S.Tm,'LineStyle',ternary(r==1,'-','--'),'Color',C(b,:));
        xc = cp_crossing(S.x,S.Tm,40); if ~isnan(xc), plot(xc*mm,40,'ko','MarkerFaceColor','y'); end
    end, end
    yline(40,'r-','40 ^oC'); yline(37,'k:','37 ^oC'); grid on; xlabel('x [mm]'); ylabel('T_m [^oC]');
    title(sprintf('Scenario "%s": T_m(x), markers = x(T_m=40 ^oC)',sc{1}));
end
legend(hh,{'T_w const, N','T_w const, CY','q_w const, N','q_w const, CY','counter, N','counter, CY','co, N','co, CY'},'Location','southeast','NumColumns',2);
sv(f,'F4_case_e_cup_mixing.png'); close(f);

f = newfig([50 50 1100 400]); tiledlayout(1,2,'TileSpacing','compact');
nexttile; plot(Qmap,dmap(:,1)*mm,'o-',Qmap,dmap(:,2)*mm,'s-', Qmap,dmap(:,5)*mm,'o--',Qmap,dmap(:,6)*mm,'s--'); grid on;
xlabel('blood flow [mL/min]'); ylabel('length [mm]');
legend('L_{max}: T_m=40 ^oC (N)','L_{max}: T_m=40 ^oC (CY)','L_{min}: T_m=37 ^oC (N)','L_{min}: T_m=37 ^oC (CY)','Location','northwest');
title('Rewarming, T_w = 42 ^oC: thermal length window');
nexttile; semilogy(Qmap,dmap(:,3),'o-',Qmap,dmap(:,4),'s-'); grid on;
xlabel('blood flow [mL/min]'); ylabel('L_{max} [m]'); legend('Newtonian','Carreau-Yasuda');
title('Pressure-limited length (p_{in} \leq 150 mmHg)');
sv(f,'F4b_case_e_f_design_map.png'); close(f);

%% Fig 5  (f) pressure
f = newfig([50 50 1100 380]); tiledlayout(1,2,'TileSpacing','compact');
nexttile; hold on; k=0;
for b=1:4, for r=1:2
    k=k+1; S = RES.(['cool_' rh{r} '_' bcs{b}]);
    plot(S.x*mm,prs(S),'LineStyle',ternary(r==1,'-','--'),'Color',C(b,:)); end, end
S = RES.full; plot(S.x*mm,prs(S),'k-','LineWidth',2);
xlabel('x [mm]'); ylabel('p [mmHg]'); grid on; title('Pressure distribution (cooling)');
legend({'T,N','T,CY','q,N','q,CY','ctr,N','ctr,CY','co,N','co,CY','full model'},'NumColumns',2,'Location','northeast');
nexttile; hold on;
for r=1:2, S = RES.(['cool_' rh{r} '_T']); plot(S.x*mm,S.tau_w,'Color',C(r,:)); end
S = RES.tdep_T; plot(S.x*mm,S.tau_w,'Color',C(3,:)); S = RES.full; plot(S.x*mm,S.tau_w,'k');
set(gca,'YScale','log'); grid on; xlabel('x [mm]'); ylabel('\tau_w [Pa]'); yline(150,'r-','150 Pa limit');
legend('Newtonian','CY','CY+\mu(T)','full (clot)','Location','northwest'); title('Wall shear stress');
sv(f,'F5_case_f_pressure.png'); close(f);

%% Fig 6  (g) clot
f = newfig([50 50 1200 700]); tiledlayout(2,3,'TileSpacing','compact');
S0 = RES.cool_carreau_yasuda_T; S1 = RES.clot_carreau_yasuda_T;
nexttile; plot(S1.x*mm,S1.R*mm,'k',S1.x*mm,-S1.R*mm,'k'); hold on; yline(base.R0*mm,'k:'); yline(-base.R0*mm,'k:');
xlabel('x [mm]'); ylabel('r [mm]'); title('Lumen with wall deposit'); grid on; ylim([-1 1]*base.R0*mm*1.1);
nexttile; plot(S0.x*mm,S0.Tm,S1.x*mm,S1.Tm); hold on;
plot(S0.x*mm,S0.Twall,'--',S1.x*mm,S1.Twall,'--'); grid on; xlabel('x [mm]'); ylabel('T [^oC]');
legend('T_m clean','T_m clot','T_{wall} clean','T_{wall} clot'); title('Temperatures (CY, const T_w)');
nexttile; plot(S0.x*mm,prs(S0),S1.x*mm,prs(S1)); grid on; xlabel('x [mm]'); ylabel('p [mmHg]');
legend('clean','clot'); title('Pressure');
nexttile; semilogy(S0.x*mm,S0.tau_w,S1.x*mm,S1.tau_w); grid on; xlabel('x [mm]'); ylabel('\tau_w [Pa]'); title('Wall shear stress');
nexttile; X = repmat(S1.x*mm,numel(S1.eta),1); Y = S1.eta*S1.R*mm;
contourf(X,Y,S1.Tfield,30,'LineColor','none'); colorbar; xlabel('x [mm]'); ylabel('r [mm]'); title('T(r,x) with clot [^oC]');
nexttile; yyaxis left; plot(dl,sev(:,2),'o-'); ylabel('\Delta p [mmHg]');
yyaxis right; plot(dl,sev(:,3),'s-'); ylabel('max \tau_w [Pa]'); xlabel('deposit severity \delta = 1-R_{min}/R_0'); grid on;
title('Severity sweep');
sv(f,'F6_case_g_clot.png'); close(f);

%% Fig 7  (h) viscous dissipation
f = newfig([50 50 1000 380]); tiledlayout(1,2,'TileSpacing','compact');
nexttile; loglog(Qlist,dis(:,1),'o-'); grid on; xlabel('blood flow [mL/min]'); ylabel('Br'); title('Brinkman number');
nexttile; loglog(Qlist,abs(dis(:,2)),'o-',Qlist,dis(:,3),'s-',Qlist,dis(:,4),'k--'); grid on;
xlabel('blood flow [mL/min]'); ylabel('\Delta T [K]');
legend('|\Delta T_m(L)| due to dissipation (const T_w)','adiabatic: T_m(L)-T_{in}','\Delta p/(\rho c_p)','Location','northwest');
title('Effect of viscous dissipation (CY)');
sv(f,'F7_case_h_dissipation.png'); close(f);

%% Fig 8  (i) temperature-dependent viscosity
f = newfig([50 50 1200 380]); tiledlayout(1,3,'TileSpacing','compact');
S0 = RES.cool_carreau_yasuda_T; S1 = RES.tdep_T;
nexttile; contourf(S1.x*mm, S1.eta*base.R0*mm, log10(S1.mufield*1e3), 30,'LineColor','none'); cb=colorbar; cb.Label.String='log_{10} \mu [mPa s]';
xlabel('x [mm]'); ylabel('r [mm]'); title('\mu(\gamma'',T) field, CY + \mu(T)');
nexttile; plot(S0.x*mm,prs(S0),S1.x*mm,prs(S1)); grid on; xlabel('x [mm]'); ylabel('p [mmHg]');
legend('CY','CY + \mu(T)'); title('Pressure');
nexttile; plot(S0.x*mm,S1.Tm-S0.Tm); hold on; plot(S0.x*mm,S1.tau_w./S0.tau_w-1); grid on;
xlabel('x [mm]'); legend('T_m difference [K]','\tau_w relative change [-]'); title('Effect of \mu(T)');
sv(f,'F8_case_i_Tdep_viscosity.png'); close(f);

%% Fig 9 full model overview
S = RES.full;
f = newfig([50 50 1200 380]); tiledlayout(1,3,'TileSpacing','compact');
nexttile; X = repmat(S.x*mm,numel(S.eta),1); Y = S.eta*S.R*mm;
contourf(X,Y,S.Tfield,30,'LineColor','none'); colorbar; xlabel('x [mm]'); ylabel('r [mm]'); title('Full model: T(r,x) [^oC]');
nexttile; plot(S.x*mm,S.Tm,S.x*mm,S.Twall,S.x*mm,S.Tc); yline(40,'r:'); yline(37,'k:'); grid on;
legend('T_m','T_{wall}','T_c (counter)'); xlabel('x [mm]'); ylabel('T [^oC]'); title('Temperatures');
nexttile; semilogx(S.HI_eta(1:end-1), S.eta(1:end-1)); grid on; xlabel('HI [%]'); ylabel('r/R');
title(sprintf('Hemolysis index per streamline (mean %.2e %%)',S.HI));
sv(f,'F9_full_model.png'); close(f);
end

function o = ternary(c,a,b), if c, o=a; else, o=b; end, end

function y = prs(S)   % gauge pressure [mmHg]
y = (S.p.p_out_mmHg*133.322 + S.dp - S.p_rel)/133.322;
end

function f = newfig(pos)
f = figure('Visible','off','Position',pos,'Color','w');
try, theme(f,'light'); catch, end
end
