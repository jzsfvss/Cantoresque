% Returns phi(adr(1))*...*phi(adr(L)).
% adr is in one/two-row vector format.

function prd = Phi(adr, phi)

sz = size(adr);
L = sz(2);
prd = 1;

if (sz(1) == 1)
	for (j = 1:L)
		prd = prd*phi(adr(j));
	end
else
	for (j = 1:L)
		prd = prd*(phi(adr(1,j))^adr(2,j));
	end
end
