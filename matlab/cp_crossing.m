function xc = cp_crossing(x, y, level)
%CP_CROSSING  First x at which y(x) crosses LEVEL (linear interpolation).
%   Returns NaN if y never reaches the level inside the domain.
d  = y - level;
s0 = sign(d(1));
if s0 == 0, xc = x(1); return; end
k = find(sign(d) ~= s0, 1, 'first');
if isempty(k), xc = NaN; return; end
xc = x(k-1) + (x(k)-x(k-1)) * d(k-1)/(d(k-1)-d(k));
end
