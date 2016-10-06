% Center and radius of the general circle for multi-map IFS (degree N=1).
% Input: 	phi, p
% Output:	[c, r]

function res = CircleGen(phi, p)

% arithmetic center = center of mass
cr2 = [ mean(p), Rho(p, mean(p)) ];

% harmonic center

n = max(size(p));
fpj = Rho(p, p(1));
for j = 2:n
	fpj = [ fpj, Rho(p, p(j)) ];
end
ifpj = 1./fpj;
hc = sum(ifpj.*p)/sum(ifpj);
cr3 = [ hc, Rho(p, hc) ];

% fminsearch
[c1, rho1] = fminsearch(@(c)max(abs(p-c)), cr2(1));
cr1 = [c1, Rho(p, c1)];

vcr = [ cr1; cr2; cr3 ];
minr = min(vcr(:,2));
for (j = 1:3)
	if (vcr(j,2) == minr)
		c = vcr(j,1);
		rho = vcr(j,2);
	end
end

lam = max(abs(phi));
mu = max(abs(1-phi));
r = (mu/(1-lam))*rho;

res = [c, r];
