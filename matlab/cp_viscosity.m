function mu = cp_viscosity(gd, aT, p)
%CP_VISCOSITY  Apparent viscosity mu(gamma_dot, T) [Pa s].
%   gd : shear rate [1/s] (array),  aT : temperature shift factor (array or scalar)
gd = abs(gd);
switch lower(p.rheology)
    case 'newtonian'
        mu = aT .* p.mu_N .* ones(size(gd));
    case 'powerlaw'
        g  = max(gd, 1e-6);                 % regularise at g -> 0
        mu = aT .* p.PL.K .* g.^(p.PL.n-1);
    case {'carreau-yasuda','cy'}
        c  = p.CY;
        mu = aT .* (c.muinf + (c.mu0-c.muinf).*(1+(c.lambda*gd).^c.a).^((c.n-1)/c.a));
    otherwise
        error('cp_viscosity:rheology','Unknown rheology "%s"', p.rheology);
end
end
