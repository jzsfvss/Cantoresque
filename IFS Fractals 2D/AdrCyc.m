% Cyclical permutations of an address.
% adr is converted to a row vector
% adrc lists row-by-row the permutations of adr

function adrc = AdrCyc(adr)

sz = size(adr);
if (sz(1) > 1) adr=AdrConv(adr); end
L = length(adr);

adrc = adr;
for k = 2:L
	adrc = [ adrc; [ adr(k:L), adr(1:(k-1)) ]];
end
