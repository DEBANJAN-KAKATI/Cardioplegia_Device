function gd = cp_shear_rate(tau, aT, p)
%CP_SHEAR_RATE  Invert the constitutive law  tau = mu(gd,T)*gd  for gd.
%   Vectorised; Carreau-Yasuda is solved with a bracketed Newton iteration
%   in log(gd) (always converges because tau(gd) is monotone for n>0).
tau = abs(tau(:));
aT  = aT(:) .* ones(size(tau));
gd  = zeros(size(tau));
pos = tau > 0;
switch lower(p.rheology)
    case 'newtonian'
        gd(pos) = tau(pos) ./ (aT(pos)*p.mu_N);
    case 'powerlaw'
        gd(pos) = (tau(pos)./(aT(pos)*p.PL.K)).^(1/p.PL.n);
    case {'carreau-yasuda','cy'}
        c  = p.CY;  t = tau(pos);  a = aT(pos);
        lo = log(t./(a*c.mu0));            % mu <= mu0   -> gd >= t/mu0
        hi = log(t./(a*c.muinf));          % mu >= muinf -> gd <= t/muinf
        z  = 0.5*(lo+hi);
        m  = (c.n-1)/c.a;
        for it = 1:100
            g  = exp(z);
            X  = (c.lambda*g).^c.a;
            B  = 1 + X;
            mu = c.muinf + (c.mu0-c.muinf)*B.^m;
            f  = log(a.*g.*mu) - log(t);              % increasing in z
            % update bracket
            hi(f>0) = z(f>0);   lo(f<=0) = z(f<=0);
            s  = 1 + (c.mu0-c.muinf)*(c.n-1)*B.^(m-1).*X ./ mu;   % dlog(tau)/dlog(g)
            zn = z - f./s;
            bad = zn<=lo | zn>=hi | ~isfinite(zn);
            zn(bad) = 0.5*(lo(bad)+hi(bad));          % bisection fallback
            if max(abs(f)) < 1e-13, break; end
            z = zn;
        end
        gd(pos) = exp(z);
    otherwise
        error('cp_shear_rate:rheology','Unknown rheology "%s"', p.rheology);
end
end
