% Converts address to a string.
% adr is in two-row format
% Detects periodicity.

function str = AdrConvStr(adr)

adr2 = AdrConv(adr);
l = length(adr2);
[ p, adr2 ] = SeqPeriod(adr2);
adr = AdrConv(adr2);
if (p < l-3) tper=1; else tper=0; end

sz = size(adr);
str = [ num2str(adr(1,1)) '^' num2str(adr(2,1)) ];
mxlen = min(sz(2),20);
for (k = 2:mxlen)
	str = [ str ' ' num2str(adr(1,k)) '^' num2str(adr(2,k)) ];
end
if (sz(2) > 20) str = [ str '...' ]; end
if (tper) str = [ str ' ' 'REP' ]; end
