% Calculates the similarity dimension of a similude IFS fractal.
% phi = IFS factors (row vector of complex numbers)

function d = DimSim(phi)

d = fsolve(@(x) sum(abs(phi).^x)-1, 1.5);
