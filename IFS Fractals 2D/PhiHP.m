% Returns phi(adr(1))*...*phi(adr(L)) (high precision).
% adr is in one/two-row vector format.

function prd = PhiHP(adr, phi)

sz = size(adr);
L = sz(2);
prd = mp('1');

if (sz(1) == 1)
	for (j = 1:L)
		prd = mp(prd*phi(adr(j)));
	end
else
	for (j = 1:L)
		prd = mp(prd*(phi(adr(1,j))^adr(2,j)));
	end
end
