% Collects the like powers in an address and creates a two-row matrix.
% OR: It creates a row vector from a matrix representation.

function adr2 = AdrConv(adr)

sz = size(adr);

if (sz(1) == 1)

A = length(adr);
adr2 = [ adr(1); 1 ];
n = 1; % mtx col counter
k = 1; % vec col counter
while (k < A)
	k = k+1;
	if (adr(k) == adr2(1,n))
		adr2(2,n) = adr2(2,n)+1;
	else
		n = n+1;
		adr2 = [ adr2, [adr(k); 1]];
	end
end

else

adr2 = adr(1,1)*ones(1,adr(2,1));
for k = 2:sz(2)
	adr2 = [ adr2, adr(1,k)*ones(1,adr(2,k)) ];
end

end
