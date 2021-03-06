% Apply an address backwards (seq of maps) to a point.
% z = 	the point
% adr =	a mtx with two rows
% 	1st row: map index
% 	2nd row: map power
% 	given in reverse the order of map application
% phi = factors of contraction (row vector of complex numbers)
% p = 	fixed points (row vector of complex numbers)
% Note: max(adr(1,:)) <= length(phi) needed.

function res = AdrB(z, adr, phi, p)

[ha, la] = size(adr);

for j = la:(-1):1
	k = adr(1,j);
	xp = adr(2,j);
	z = p(k) + (phi(k)^xp)*(z-p(k));
end

res = z;
