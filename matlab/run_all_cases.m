function T = run_all_cases()
%RUN_ALL_CASES  Solves problem statement 4(a)-(i) and writes figures + tables
%               to ../results.
%
%   Scenario "cool"  : locally heated blood (T_in = 42 C) cooled by the HX.
%                      -> x where Tm <= 40 C = MIN length for safety,
%                         x where Tm <= 37 C = length to reach normothermia.
%   Scenario "rewarm": blood (32 C) rewarmed by 42 C water.
%                      -> x where Tm >= 40 C = MAX allowable length (4e).
close all;
here   = fileparts(mfilename('fullpath'));
outdir = fullfile(here,'..','results');
if ~exist(outdir,'dir'), mkdir(outdir); end
diary(fullfile(outdir,'run_log.txt')); diary on;
t0 = tic;
base = cp_params();
scen.cool   = struct('T_in',42,'T_wall',32,'q_wall', 3000,'Tc_in',28);
scen.rewarm = struct('T_in',32,'T_wall',42,'q_wall',-5000,'Tc_in',42);
rheos = {'newtonian','carreau-yasuda'};
bcs   = {'T','q','coupled-counter','coupled-co'};
rows  = {};  RES = struct();

%% ---------------- (a)-(e): 2 scenarios x 2 rheologies x 4 BCs ----------
for sn = fieldnames(scen)'
    sc = scen.(sn{1});
    for r = 1:numel(rheos)
        for b = 1:numel(bcs)
            p = setscen(base, sc); p.rheology = rheos{r};
            p = setbc(p, bcs{b});
            S = cp_solve(p);
            key = matlab.lang.makeValidName([sn{1} '_' rheos{r} '_' bcs{b}]);
            RES.(key) = S;
            rows(end+1,:) = metrics(S, sn{1}, sprintf('%s | %s', bcs{b}, rheos{r})); %#ok<AGROW>
            fprintf('%-8s %-16s %-16s Tm(L)=%6.3f  dp=%5.2f mmHg  [%.1fs]\n', ...
                sn{1}, rheos{r}, bcs{b}, S.Tm(end), S.dp_mmHg, toc(t0));
        end
    end
end

%% ---------------- (c) extra: power-law comparison ----------------------
p = setscen(base, scen.cool); p.rheology = 'powerlaw'; p = setbc(p,'T');
RES.cool_powerlaw_T = cp_solve(p);
rows(end+1,:) = metrics(RES.cool_powerlaw_T,'cool','T | powerlaw');

%% ---------------- (d) coolant-flow sweep --------------------------------
Qc_list = [0.25 0.5 1 2 5 10];          % L/min
sweep = zeros(numel(Qc_list),4);
for k = 1:numel(Qc_list)
    for a = 1:2
        p = setscen(base, scen.cool); p.rheology = 'carreau-yasuda';
        p = setbc(p, bcs{2+a}); p.Qc_total = Qc_list(k)/60/1000;
        S = cp_solve(p);
        sweep(k,a)   = S.Tm(end);
        sweep(k,a+2) = S.Tc(1);         % counter: coolant outlet; co: inlet
    end
end
fprintf('(d) coolant sweep done [%.1fs]\n',toc(t0));

%% ---------------- (g) wall deposit / clot -------------------------------
for b = 1:numel(bcs)
    for r = 1:numel(rheos)
        p = setscen(base, scen.cool); p.rheology = rheos{r};
        p = setbc(p, bcs{b}); p.clot.on = true;
        S = cp_solve(p);
        key = matlab.lang.makeValidName(['clot_' rheos{r} '_' bcs{b}]);
        RES.(key) = S;
        rows(end+1,:) = metrics(S,'cool+clot',sprintf('%s | %s', bcs{b}, rheos{r})); %#ok<AGROW>
    end
end
% severity sweep (CY, constant wall temperature)
dl = [0 0.2 0.4 0.5 0.6 0.7];  sev = zeros(numel(dl),4);
for k = 1:numel(dl)
    p = setscen(base, scen.cool); p.rheology='carreau-yasuda'; p = setbc(p,'T');
    p.clot.on = dl(k)>0; p.clot.delta = dl(k);
    S = cp_solve(p);
    sev(k,:) = [S.Tm(end) S.dp_mmHg max(S.tau_w) S.HI];
end
fprintf('(g) clot cases done [%.1fs]\n',toc(t0));

%% ---------------- (h) viscous dissipation -------------------------------
Qlist = [100 200 300 500 1000 2000];      % mL/min total
dis = zeros(numel(Qlist),5);
for k = 1:numel(Qlist)
    p = setscen(base, scen.cool); p.rheology='carreau-yasuda'; p = setbc(p,'T');
    p.Q_total = Qlist(k)*1e-6/60;
    S0 = cp_solve(p);  p.viscDiss = true;  S1 = cp_solve(p);
    pa = p; pa.bc='q'; pa.q_wall=0; pa.T_in=37; Sa = cp_solve(pa);   % adiabatic
    U  = S0.Qt/(pi*p.R0^2);
    Br = p.mu_N*U^2/(p.k*abs(p.T_in-p.T_wall));
    dis(k,:) = [Br, S1.Tm(end)-S0.Tm(end), Sa.Tm(end)-37, Sa.dp/(p.rho*p.cp), max(S1.Tmax-S0.Tmax)];
end
p = setscen(base, scen.cool); p.rheology='carreau-yasuda'; p=setbc(p,'T'); p.viscDiss=true;
RES.diss_cy_T = cp_solve(p);
rows(end+1,:) = metrics(RES.diss_cy_T,'cool','T | CY + dissipation');
fprintf('(h) dissipation done [%.1fs]\n',toc(t0));

%% ---------------- (i) temperature-dependent viscosity -------------------
for b = 1:numel(bcs)
    p = setscen(base, scen.cool); p.rheology='carreau-yasuda'; p = setbc(p,bcs{b});
    p.Tdep = true;
    S = cp_solve(p);
    key = matlab.lang.makeValidName(['tdep_' bcs{b}]);
    RES.(key) = S;
    rows(end+1,:) = metrics(S,'cool',sprintf('%s | CY + mu(T)', bcs{b})); %#ok<AGROW>
end
% full model: CY + mu(T) + dissipation + clot + counter-current coupling
p = setscen(base, scen.cool); p.rheology='carreau-yasuda'; p = setbc(p,'coupled-counter');
p.Tdep = true; p.viscDiss = true; p.clot.on = true;
RES.full = cp_solve(p);
rows(end+1,:) = metrics(RES.full,'cool+clot','coupled-counter | CY+mu(T)+diss (FULL)');
fprintf('(i) T-dependent viscosity done [%.1fs]\n',toc(t0));

%% ---------------- design map (4e + 4f): limits vs blood flow -----------
Qmap = [100 150 200 300 400 500 600];
dmap = nan(numel(Qmap),6);
for k = 1:numel(Qmap)
    for r = 1:2
        p = setscen(base, scen.rewarm); p.rheology = rheos{r}; p = setbc(p,'T');
        p.Q_total = Qmap(k)*1e-6/60; p.L = 0.6; p.Nx = 600;
        S = cp_solve(p);
        dmap(k,r)   = cp_crossing(S.x,S.Tm,p.T_safe);                   % L_max thermal
        dpdx_mean   = S.dp/p.L;
        dmap(k,r+2) = (p.p_max_mmHg - p.p_out_mmHg)*133.322/dpdx_mean;  % L_max pressure
        dmap(k,r+4) = cp_crossing(S.x,S.Tm,p.T_target);                 % L_min target
    end
end
fprintf('design map done [%.1fs]\n',toc(t0));

%% ---------------- tables --------------------------------------------------
hdr = {'scenario','case','Tm_out_C','x_Tm40_mm','x_Tm37_mm','Tmax_local_C','Tmin_local_C', ...
       'dp_mmHg','p_in_mmHg','p_ok','tau_w_max_Pa','tau_ok','HI_percent','Qheat_W','Nu_out','energy_err'};
T = cell2table(rows,'VariableNames',hdr);
writetable(T, fullfile(outdir,'summary_cases.csv'));
disp(T);
writematrix([Qc_list' sweep], fullfile(outdir,'d_coolant_sweep.csv'));
writematrix([dl' sev], fullfile(outdir,'g_clot_severity.csv'));
writematrix([Qlist' dis], fullfile(outdir,'h_dissipation.csv'));
writematrix([Qmap' dmap], fullfile(outdir,'e_design_map.csv'));
save(fullfile(outdir,'all_results.mat'),'RES','T','sweep','Qc_list','sev','dl','dis','Qlist','dmap','Qmap','-v7.3');

%% ---------------- figures ------------------------------------------------
make_figures(RES, outdir, base, sweep, Qc_list, sev, dl, dis, Qlist, dmap, Qmap);
fprintf('\nTotal time %.1f s. Results in %s\n', toc(t0), outdir);
diary off;
end

% =========================================================================
function p = setscen(p, sc)
f = fieldnames(sc); for i=1:numel(f), p.(f{i}) = sc.(f{i}); end
end

function p = setbc(p, bc)
switch bc
    case 'T', p.bc = 'T';
    case 'q', p.bc = 'q';
    case 'coupled-counter', p.bc = 'coupled'; p.flow_arrangement = 'counter';
    case 'coupled-co',      p.bc = 'coupled'; p.flow_arrangement = 'co';
end
end

function row = metrics(S, scen, label)
p = S.p;
x40 = cp_crossing(S.x,S.Tm,p.T_safe)*1e3;
x37 = cp_crossing(S.x,S.Tm,p.T_target)*1e3;
row = {scen, label, S.Tm(end), x40, x37, max(S.Tmax), min(S.Tmin), S.dp_mmHg, ...
       S.p_in_mmHg, S.p_in_mmHg <= p.p_max_mmHg, max(S.tau_w), max(S.tau_w) <= p.tau_max, ...
       S.HI, S.Qheat_total, S.Nu(end), S.energy_err};
end
