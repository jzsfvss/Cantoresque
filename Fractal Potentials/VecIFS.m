% Returns the vector field IFS iterates to a certain level N.
% Z =	grid points
% phi = factors of contraction (row vector of complex numbers)
% p = 	fixed points (row vector of complex numbers)
% prb =	probabilities (row vector of real numbers in (0,1))
% N = 	iteration level (natural number)
% typ = type of the initial flow (superposition of eddies, or the character eddy)

function res = VecIFS(Z, phi, p, prb, N, typ)

n = length(phi);
VecZ = 0;

if (N == 1)

% Tranfer from initial level:
if (typ == 1)
	for (k = 1:n)
		Zk = p(k) + (1/phi(k))*(Z - p(k));
		VecZk = VecVor(Zk, phi(k), p(k));
		VecZ = VecZ + prb(k)*(phi(k)/(abs(phi(k))^2))*VecZk;
	end
else
	pave = sum(p)/n;
	char = sum(prb.*log(phi));
	VecZ = char*(Z-pave)./(abs(Z-pave).^2);
end

else

% Transfer to next level:
for (k = 1:n)
	Zk = p(k) + (1/phi(k))*(Z - p(k));
	VecZk = VecIFS(Zk, phi, p, prb, N-1, typ);
	VecZ = VecZ + prb(k)*(phi(k)/(abs(phi(k))^2))*VecZk;
end

end

res = VecZ;
