% Returns the stream IFS iterates of S to a certain level N.
% Z =	data points
% phi = factors of contraction (row vector of complex numbers)
% p = 	fixed points (row vector of complex numbers)
% prb =	probabilities (row vector of real numbers in (0,1))
% N = 	iteration level (natural number)

function res = PsiIFS(Z, phi, p, prb, N)

n = length(phi);
% [lr,lc] = size(Z);
% PsiZ = zeros(lr,lc);
PsiZ = 0;
% ref = 0.5;
ref = p + 1;

if (N == 1)

% Tranfer from initial level:
for (k = 1:n)
	Zk = p(k) + (1/phi(k))*(Z - p(k));
%	PsiZk = angle((Zk-p(k))/(ref-p(k)));
	PsiZk = PsiVor(Zk, phi(k), p(k), ref(k));
	PsiZ = PsiZ + prb(k)*PsiZk;
end

else

% Transfer to next level:
for (k = 1:n)
	Zk = p(k) + (1/phi(k))*(Z - p(k));
	PsiZk = PsiIFS(Zk, phi, p, prb, N-1);
	PsiZ = PsiZ + prb(k)*PsiZk;
end

end

res = PsiZ;
