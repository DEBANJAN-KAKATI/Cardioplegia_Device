# Cardioplegia Heat Exchanger: 2-D Laminar Thermo-Fluid Model (MATLAB)

Complete solution of problem statement 4 (a)–(i): steady, laminar, axisymmetric
heat transfer of blood in the tubes of a cardioplegia heat exchanger. The model uses the
lubrication approximation, with Newtonian, Carreau–Yasuda and power-law rheology, three
wall conditions (constant T, constant q, coupled coolant circuit), a cup-mixing safety
length, the pressure distribution, wall deposits (clots), viscous dissipation and
temperature-dependent viscosity. Everything is verified against analytical solutions.

> **Before you publish:** the geometry and flow values in `matlab/cp_params.m` are
> *generic placeholders* (40 tubes, ID 1.5 mm, L = 150 mm, 300 mL/min). The slides contain
> no dimensions. Put in the real device's tube count, ID, length, wall material, water flow
> and h_c, then re-run. Everything else adapts automatically.

---

## 1. How to run

```matlab
cd matlab
run_everything          % steady + transient validation and all cases (~20 min); writes ../results/*
```

| File | Purpose |
|---|---|
| `cp_params.m` | **All** parameters (geometry, blood/water properties, rheology, BCs, limits, grid) |
| `cp_solve.m` | Core 2-D marching solver (energy + lubrication momentum + coolant circuit) |
| `cp_velocity.m` | Local fully-developed velocity for any rheology at a given flow rate |
| `cp_shear_rate.m`, `cp_viscosity.m` | Constitutive laws (Newtonian / Carreau–Yasuda / power law, × a_T(T)) |
| `cp_radius.m` | Lumen radius R(x) with an optional wall deposit (clot) |
| `cp_crossing.m` | x where T_m crosses a level (safety lengths, item e) |
| `cp_validate.m` | Verification suite (V1–V6) |
| `run_all_cases.m`, `make_figures.m` | Cases (a)–(i), sweeps, design maps, figures, tables |
| `cp_solve_transient.m` | Transient blood + wall + water model with time-varying inputs (§5) |
| `cp_validate_transient.m`, `run_transient_cases.m` | Transient verification (TV1–TV4) and digital-twin scenarios (D1–D5) |

Single run example:

```matlab
p = cp_params;  p.rheology = 'carreau-yasuda';  p.bc = 'coupled';
p.flow_arrangement = 'counter';  p.Tdep = true;  p.viscDiss = true;
S = cp_solve(p);
plot(S.x, S.Tm)                         % cup-mixing temperature
cp_crossing(S.x, S.Tm, 40)              % x where Tm = 40 C
```

Outputs in `results/`: `summary_cases.csv` (all cases), `d_coolant_sweep.csv`,
`e_design_map.csv`, `g_clot_severity.csv`, `h_dissipation.csv`, `all_results.mat` (full
fields), figures `F1…F9`, validation figures `V2`, `V6`, and `run_log.txt`.
(The `.mat` files are ~20 MB each and are not tracked in git; `run_everything` regenerates them.)

---

## 2. Mathematical formulation

### 2.1 Assumptions
Steady flow. Laminar (Re ≈ 30). Axisymmetric tube of radius R(x) with R'(x) ≪ 1, so the
**lubrication approximation** holds. Hydrodynamically developed: the entrance length is
0.05 Re D ≈ 2 mm, which is less than 2 % of L. Thermally developing (Graetz). Axial
conduction is neglected because Pe = Re·Pr ≈ 800. ρ, c_p and k are constant; μ = μ(γ̇, T).

### 2.2 Momentum (lubrication approximation)
$$\frac{1}{r}\frac{\partial}{\partial r}\left(r\,\tau_{rx}\right)=\frac{dp}{dx}
\;\Rightarrow\; \tau(r,x)=\tau_w(x)\,\frac{r}{R(x)},\qquad
\tau=\mu(\dot\gamma,T)\,\dot\gamma,\quad \dot\gamma=-\frac{\partial u}{\partial r}$$
$$u(r,x)=\int_r^{R}\dot\gamma\,dr',\qquad 2\pi\!\int_0^R u\,r\,dr = Q_t
\;\;\Rightarrow\;\;\tau_w(x),\quad \frac{dp}{dx}=-\frac{2\tau_w}{R}$$
The radial velocity v follows from continuity and is carried by the stream function (§2.4).

### 2.3 Rheology
* Newtonian: $\mu=\mu_N\,a_T$
* Carreau–Yasuda (Boyd 2007 / Cho & Kensey 1991): $\mu=a_T\left[\mu_\infty+(\mu_0-\mu_\infty)\left(1+(\lambda\dot\gamma)^a\right)^{(n-1)/a}\right]$, with μ₀ = 0.056, μ∞ = 0.00345 Pa s, λ = 1.902 s, n = 0.22, a = 1.25
* Power law (comparison only): $\mu=a_T K\dot\gamma^{\,n-1}$
* Temperature dependence (item i): $a_T=\exp[-\beta(T-T_{ref})]$, β = 0.021 K⁻¹ (≈2 %/K)

### 2.4 Energy equation (items a–i)
$$\rho c_p\left(u\frac{\partial T}{\partial x}+v\frac{\partial T}{\partial r}\right)=\frac{1}{r}\frac{\partial}{\partial r}\left(k r\frac{\partial T}{\partial r}\right)+\underbrace{\mu\dot\gamma^{2}}_{\text{(h) dissipation}}$$
In the mapped coordinate η = r/R(x), which follows a clotted wall exactly, this becomes
$$\rho c_p\left(u\,T_x-\frac{G_x}{R^2\eta}T_\eta\right)=\frac{1}{R^2\eta}\left(k\eta T_\eta\right)_\eta+\mu\dot\gamma^2,
\qquad G(\eta,x)=R^2\!\int_0^\eta u\,\eta'\,d\eta'.$$
Here G is the Stokes stream function divided by 2π. The G_x term is the exact form of the
v-convection: it vanishes for self-similar profiles and is non-zero under a clot or when
μ(T) changes.

**Inlet / symmetry:** T(r,0) = T_in; ∂T/∂r = 0 at r = 0.

**Wall conditions**

| Case | Condition at the blood–wall interface |
|---|---|
| (a) constant wall temperature | $T(R_0,x)=T_w$ |
| (b) constant heat flux | $-k\,\partial T/\partial r=q_w$ (q_w > 0 cools the blood) |
| (d) coupled coolant | $-k\,\partial T/\partial r = U\,(T_{wall}-T_c(x))$ |

For (d), the second energy equation is the coolant circuit:
$$\dot m_c c_{p,c}\frac{dT_c}{dx}=\pm\,2\pi R\,U\,(T_{wall}-T_c)$$
The + sign is co-current; the − sign is counter-current, with T_c(L) = T_{c,in}, solved by
secant shooting. The series resistance is
$\frac1U=\frac{R\ln(R_0/R)}{k_{clot}}+\frac{R\ln(R_o/R_0)}{k_{wall}}+\frac{R}{R_o h_c}$.
The first term is the thrombus layer (item g). In (a) the clot also turns the Dirichlet
condition into a Robin condition.

### 2.5 Derived quantities
* **Cup-mixing temperature (e):** $T_m(x)=\dfrac{\int_0^R uTr\,dr}{\int_0^R u r\,dr}$. The safety
  length is x* where T_m = 40 °C (and T_m = 37 °C). The local maximum T is also checked against
  42 °C (protein damage).
* **Nusselt number:** $Nu=2Rq_w/[k(T_m-T_{wall})]$
* **Pressure (f):** $p(x)=p_{out}+\int_x^L 2\tau_w/R\,dx'$, with constraints p_in ≤ 150 mmHg and τ_w ≤ 150 Pa
* **Hemolysis:** Giersiepen power law $HI=C\,t^{\alpha}\tau^{\beta}$, accumulated along every
  streamline with the linearised (Garon–Fariñas) form, then flow-averaged.
* **Global energy balance** per run: $\rho c_pQ\Delta T_m=-\!\int2\pi Rq_w\,dx+\!\iint\mu\dot\gamma^2dA\,dx$.
  The reported error is below 2·10⁻⁴ in every case.

### 2.6 Numerical method
* Finite volume in η on a wall-clustered grid (η = sin(πs/2), 101 nodes)
* Fully implicit marching in x on an inlet-clustered grid (x = L s², 500 steps), which
  resolves the Graetz entrance singularity
* Sparse direct solve at each step; the coolant temperature is solved implicitly with the
  blood temperature
* Picard iteration when μ = μ(T)
* Carreau–Yasuda: τ = μ(γ̇)γ̇ is inverted with bracketed Newton iteration in log γ̇, and τ_w
  comes from a monotone scalar root of Q(τ_w) = Q_t
* Run time on a laptop: under 1 s (Newtonian or CY), up to ~1–2 min (CY + μ(T) +
  counter-current shooting). This is fast enough for optimisation and digital-twin use.

---

## 3. Verification (`cp_validate.m`)

| Test | Result |
|---|---|
| V1 Poiseuille / power-law velocity and τ_w vs analytical | error ≤ 4·10⁻⁵ |
| V1 CY with μ₀ = μ∞ reduces to Poiseuille | 3·10⁻⁵ |
| V2 Graetz problem (independent eigen-series, 80 modes): θ_m(x) | max error 5.5·10⁻⁴ |
| V2 Nu∞ (constant T_w) | **3.6567** (exact 3.6568) |
| V3 constant flux: T_m linear; Nu∞ | **4.3635** (exact 48/11 = 4.3636) |
| V4 power law, constant flux, Nu∞ = 8(5n+1)(3n+1)/(31n²+12n+1) | n = 0.5: 4.7456 vs 4.7458; n = 0.708: 4.5320 vs 4.5322 |
| V5 adiabatic wall + dissipation: ΔT_m = Δp/(ρc_p) | relative error 1.6·10⁻⁴ |
| V6 grid convergence | 1st order (as designed); default grid error 0.25 % in Nu, 3·10⁻⁴ in θ_m |

---

## 4. Results with the default (placeholder) device

The "cool" scenario is note 1: locally heated blood at 42 °C is cooled. The wall is at 32 °C
in (a), q_w = 3000 W/m² in (b), and coolant enters at 28 °C, 2 L/min in (d). The "rewarm"
scenario has blood at 32 °C and water/wall at 42 °C; this is where note (e)'s "upper limit
of length" appears.

| Item | Key finding |
|---|---|
| (a) constant T_w, Newtonian | T_m(L) = 33.35 °C; T_m ≤ 40 °C after **8.5 mm**, 37 °C at **41.8 mm**; Nu → 3.657 |
| (b) constant q_w, Newtonian | T_m falls linearly: 40 °C at 67.5 mm, outlet 37.55 °C; Nu → 4.367 |
| (c) Carreau–Yasuda | Blunter profile. τ_w and Δp **+9.6 %**, Nu∞ +1 % (3.691 / 4.409), T_m changes < 0.03 K. A power law *under*-predicts τ_w by 16 % relative to CY. |
| (d) coupled coolant | Counter-current outlet 32.53 °C vs co-current 32.91 °C at 2 L/min. At 0.25 L/min the gap is 1.5 K (34.4 vs 36.0 °C). T_wall(x) is no longer constant (fig. F3). |
| (e) safety length (rewarm, T_w = 42 °C) | **L_max = 117 mm** (N) / 116 mm (CY) at 300 mL/min. L_max grows linearly with flow: 39 mm at 100 mL/min, 235 mm at 600 mL/min (F4b). Constant-flux rewarming keeps T_m < 40 °C, yet the **local wall-blood temperature reaches 42.8 °C**, above the protein limit. |
| (f) pressure | Δp ≈ 4 mmHg and p_in ≈ 84 mmHg, well below 150 mmHg. The pressure-limited length is 1.3–8 m, so pressure does not bind for a clean device. τ_w ≈ 1.4 Pa, far below 150 Pa. |
| (g) wall deposit, δ = 40 % | Δp 4.3 → 6.5 mmHg; max τ_w 1.4 → 6.2 Pa; HI ×2.3; the thrombus insulates the wall, so T_m(L) rises 0.14 K. At δ = 70 %: Δp 26 mmHg, τ_w 48 Pa. |
| (h) viscous dissipation | Br ≈ 4·10⁻⁶; ΔT_m ≈ 3·10⁻⁵ K (adiabatic heating = Δp/ρc_p ≈ 1.5·10⁻⁴ K). Negligible at clinical flows (a useful negative result). |
| (i) μ(T) | Cooling raises viscosity, so Δp **+6.9 %**, outlet τ_w +10 %, and the 37 °C length shifts 41.4 → 42.5 mm |
| Full model (CY + μ(T) + dissipation + clot + counter-current) | T_m(L) = 32.79 °C, p_in = 86.8 mmHg, max τ_w = 6.6 Pa, HI = 1.7·10⁻⁴ % |

Figures: `F1` (a, b), `F2` (c), `F3` (d), `F4`/`F4b` (e), `F5` (f), `F6` (g), `F7` (h), `F8` (i), `F9` full model.

---

## 5. Transient model (digital-twin core)

`cp_solve_transient(p, tr)` extends the model in time. Every operating input can be any
function of time:

```matlab
p = cp_params; p.rheology = 'carreau-yasuda'; p.bc = 'coupled'; p.flow_arrangement = 'counter';
tr.t_end = 40; tr.dt = 0.05;
tr.Tc_in   = @(t) 28 - 8*(t >= 2);        % water inlet step
tr.Q_total = @(t) 5e-6*(1 + 0.4*sin(2*pi*t));   % pulsatile blood flow
S = cp_solve_transient(p, tr);
plot(S.t, S.Tm_out)
```

Time-varying inputs: `T_in(t)`, `Q_total(t)`, `Tc_in(t)`, `Qc_total(t)`, `T_wall(t)`,
`q_wall(t)`. Outputs vs time: `Tm_out`, `Tc_out`, `Tmax`, `Twall_max`, `p_in_mmHg`,
`dp_mmHg`, `Qheat_W`, `tau_w_max`, plus `Tm_xt` (T_m(x,t)), `Tc_xt` and full T(r,x)
snapshots at `tr.save_times`.

### 5.1 Equations
$$\rho c_p\left(\frac{\partial T}{\partial t}+u\,T_x-\frac{G_x}{R^2\eta}T_\eta\right)=\frac{1}{R^2\eta}(k\eta T_\eta)_\eta+\mu\dot\gamma^2\qquad\text{(blood)}$$
$$\rho_w c_w A_w\frac{\partial T_s}{\partial t}=2\pi R\left[U_1(T_{wall}-T_s)-U_2(T_s-T_c)\right]\qquad\text{(tube wall, lumped radially)}$$
$$\rho_c c_{p,c}A_c\frac{\partial T_c}{\partial t}\pm\dot m_c c_{p,c}\frac{\partial T_c}{\partial x}=2\pi R\,U_2(T_s-T_c)\qquad\text{(water)}$$
Here $1/U_1 = R_{clot}+R\ln(R_m/R_0)/k_w$ and $1/U_2 = R\ln(R_o/R_m)/k_w + R/(R_o h_c)$, with
R_m the mid-wall radius. At steady state these reduce exactly to the series resistance U of §2.4.

**Quasi-steady flow.** The momentum equation stays quasi-steady: the viscous diffusion time
R²/ν ≈ 0.2 s is much shorter than the thermal times (residence 2.1 s, radial diffusion
4.3 s). The Womersley number at 1 Hz is ≈ 1, so for pulsatile flow this is a first-order
approximation.

**Numerics.**
* Space uses the same finite-volume / upwind discretisation as the steady solver.
* Time uses BDF2, fully implicit (backward Euler optional). The sparse LU factorisation is
  reused while the flow rates are unchanged.
* Counter-current water flow is solved simultaneously with the blood, so no shooting is needed.
* With μ(T), the viscosity is lagged by one time step.
* Run time: 40 s of simulated time in about 5–10 s (61 × 200 grid).

### 5.2 Verification (`cp_validate_transient.m`)
| Test | Result |
|---|---|
| TV1 steady limit vs `cp_solve` (constant T, constant q, co-current) | 10⁻¹³–10⁻¹² K (identical discretisation) |
| TV1 counter-current vs fine reference | first-order convergence (0.034 → 0.0036 K, Nx 100 → 800) |
| TV2 discrete energy conservation (blood + wall + water) | relative residual 8·10⁻¹³ |
| TV3 Poiseuille residence-time distribution, F = 1 − (τ/2t)² | L1 error 7.4·10⁻³ → 1.1·10⁻³ (first order) |
| TV4 temporal accuracy (BDF2) | observed order ≥ 2 (≈3 at the test point) |

### 5.3 Digital-twin scenarios (`run_transient_cases.m`)
All scenarios use cooling with CY blood, counter-current water and wall heat capacity.

| Scenario | Result |
|---|---|
| D1 water inlet step 28 → 20 °C | T_m,out 32.54 → 27.14 °C. Fitted first-order-plus-dead-time (FOPDT) model: K = 0.675, dead time θ = 0.27 s, τ = 1.11 s; 95 % settled after 3.1 s |
| D2 blood inlet step +1 K | K = 0.325, θ = 1.13 s, τ = 0.72 s. These values give the transfer functions G(s) = K e^{−θs}/(τs+1) for controller design. |
| D3 blood flow +50 % | T_m,out 32.54 → 34.51 °C (τ ≈ 0.6 s); p_in 84.3 → 86.3 mmHg |
| D4 **pump stop during rewarming** (water 42 °C, flow down to 5 % for 30 s) | Outlet passes 40 °C **0.7 s after the stop**, stays above 40 °C for 31 s, peaks at 42 °C. After restart a **hot bolus** reaches the patient: **16 mL of blood above 40 °C** in total. |
| D5 pulsatile roller pump (1 Hz, ±40 %) | T_m,out oscillates ±0.22 °C, p_in ±1.6 mmHg |

Figures: `FT1_step_responses`, `FT2_flow_step`, `FT3_pump_stop_rewarming`, `FT4_pulsatile`,
plus validation figures `VT3_residence_time` and `VT4_time_convergence`.

Note: a sudden step input causes a small BDF2 overshoot in the water outlet temperature
(≈0.4 K, under 1 s). Use `tr.scheme = 'BE'` or a smoothed input if that matters.

---

## 6. Limitations and next steps
1. **Real geometry:** replace the placeholders in `cp_params.m`. If the blood side is a pleated
   or annular channel rather than tubes, only `cp_velocity` (the cross-section) needs a variant.
   The slit or annulus version is the same 1-D τ-integration.
2. h_c on the water side is a constant input. Replace it with a shell-side correlation or
   manufacturer data.
3. **Digital twin, next step:** the transient model (§5) is the physics core. What remains is
   to connect it to measured inlet/outlet temperatures and pressures (state estimation with
   an extended or ensemble Kalman filter), and to fit h_c and the clot severity from bench
   data. The FOPDT transfer functions (§5.3) can drive a model-predictive controller for
   the water temperature and flow.
4. Oxygenator (two-phase injection) is a separate problem; it is not covered here.

### Suggested publications
* **Paper 1 (modelling):** "Thermally developing laminar flow of shear-thinning,
  temperature-sensitive blood in a cardioplegia heat exchanger with conjugate coolant
  coupling". Covers items a–d, h, i, the Graetz/power-law validation, and Nu∞ for CY blood.
* **Paper 2 (application/safety):** "Safety-constrained design and operating maps for
  cardioplegia heat exchangers: temperature, pressure, wall shear and clot deposition".
  Covers items e–g, design maps, a hemolysis estimate, and real-time optimisation as a
  digital-twin core, validated against company bench data.
