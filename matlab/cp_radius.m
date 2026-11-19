function R = cp_radius(x, p)
%CP_RADIUS  Local lumen radius R(x); a Gaussian wall deposit (clot) if p.clot.on.
R = p.R0*ones(size(x));
if isfield(p,'clot') && p.clot.on
    xc = p.clot.xc*p.L;
    R  = p.R0*(1 - p.clot.delta*exp(-((x-xc)/p.clot.w).^2));
end
end
