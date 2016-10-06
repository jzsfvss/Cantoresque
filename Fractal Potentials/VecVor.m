% Calculates the velocity vector field of a vortex at certain points.
% Input: 	z (several complex numbers)
%		phi, p (a single complex number)
% Output:	complex velocity vectors at the points in z
% C0 = 2*pi

function v = VecVor(z, phi, p)

v = log(phi)*(z-p)./(abs(z-p).^2);
