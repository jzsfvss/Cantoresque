% Converts address to a string from 1/2-row format.
% Does not detect periodicity.

function str = AdrConvStr(adr)

sz = size(adr);
if (sz(1) == 1)
	adr = AdrConv(adr);
end

str = [ num2str(adr(1,1)) ];
if (adr(2,1) > 1)
		str = [ str '^' num2str(adr(2,1)) ];
end
for (k = 2:sz(2))
	str = [ str ' ' num2str(adr(1,k)) ];
	if (adr(2,k) > 1)
		str = [ str '^' num2str(adr(2,k)) ];
	end
end
