% Calculates the conjugate of a vector of IFS eddy phis.

function phi2 = EdConj(phi)

phi2 = phi.^(i*sign(angle(phi)));
