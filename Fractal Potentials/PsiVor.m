% Calculates the stream function of a vortex.
% Input: 	z, phi, p (fixed pt), rp (ref pt)
% Output:	psi value

function res = PsiVor(z, phi, p, rp)

pitch = angle(phi)/(log(abs(phi)));
vec = (z-p)/(rp-p);

res = mod(angle(vec)-pitch*log(abs(vec)), 2*pi);
