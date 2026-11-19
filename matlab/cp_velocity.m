function V = cp_velocity(R, Qt, eta, T, p)
%CP_VELOCITY  Local fully-developed (lubrication) velocity profile.
%
%   Momentum (lubrication approximation, |dR/dx| << 1):
%       (1/r) d/dr ( r tau_rx ) = dp/dx   ->   tau(r) = tau_w * r/R
%   with tau = mu(gd,T) gd ,  gd = -du/dr ,  u(R) = 0,
%   and tau_w chosen so that  2*pi*int_0^R u r dr = Qt  (flow per tube).
%
%   Inputs : R radius at this x, Qt flow per tube, eta = r/R nodes (col),
%            T local temperature profile on eta (needed if p.Tdep).
%   Output : struct with u, gd, mu, tau, Phi (viscous dissipation),
%            tau_w, dpdx.
eta = eta(:);
if p.Tdep
    aT = exp(-p.beta_mu*(T(:) - p.T_ref));
else
    aT = ones(size(eta));
end

    function [Q,u,gd] = flowrate(tw)
        gd = cp_shear_rate(tw*eta, aT, p);
        I  = cumtrapz(eta, gd);
        u  = R*(I(end) - I);                    % u = R * int_eta^1 gd d(eta)
        Q  = 2*pi*R^2*trapz(eta, u.*eta);
    end

switch lower(p.rheology)
    case 'newtonian'          % Q linear in tau_w
        Q1 = flowrate(1);  tw = Qt/Q1;
    case 'powerlaw'           % Q ~ tau_w^(1/n)
        Q1 = flowrate(1);  tw = (Qt/Q1)^p.PL.n;
    otherwise                 % Carreau-Yasuda: monotone scalar root
        tlo = 4*p.CY.muinf*min(aT)*Qt/(pi*R^3);   % Poiseuille with mu_inf
        thi = 4*p.CY.mu0  *max(aT)*Qt/(pi*R^3);   % Poiseuille with mu0
        fun = @(lt) log(flowrate(exp(lt))/Qt);
        opt = optimset('TolX',1e-13,'Display','off');
        tw  = exp(fzero(fun, [log(tlo*0.999) log(thi*1.001)], opt));
end
[~,u,gd] = flowrate(tw);
V.u     = u;
V.gd    = gd;
V.mu    = cp_viscosity(gd, aT, p);
V.tau   = tw*eta;
V.Phi   = V.tau.*gd;                 % viscous dissipation  mu*gd^2  [W/m^3]
V.tau_w = tw;
V.dpdx  = -2*tw/R;                   % [Pa/m]
V.aT    = aT;
end
