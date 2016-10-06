% Determines the period of a sequence.
% p = 	period
% ap =	periodic part of vector a

function [ p, ap ] = SeqPeriod(a)

N = length(a);
p = 0;
tst1 = 1;

while (tst1)
	p = p+1;
	tst2 = 0;
	if (a(p+1) == a(1))
		tst2 = 1;
		k = 1;
		while (tst2 && (p+k+1 <= N))
			k = k+1;
			if (a(p+k) ~= a(k)) tst2=0; end
		end
	end
	if (tst2 || (p+2 > N)) tst1=0; end
end

ap = a(1:p);
