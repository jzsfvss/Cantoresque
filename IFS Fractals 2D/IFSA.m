% Returns the IFS iterates of S to a certain level N along with the addresses wrt p(1).
% N = 	iteration level (natural number)
% A = 	max length of address
% phi = factors of contraction (row vector of complex numbers)
% p = 	fixed points (row vector of complex numbers)
% adr = addresses (|F|xA matrix of integers)

function [f, adr] = IFSA(phi, p, N, A)

n=length(phi);

if (N == 0)

f = p(1);
adr = 0;

else

% Previous level points:
[f, adr] = IFSA(phi, p, N-1, A-1);
lf = length(f);

f2 = p(1) + phi(1)*(f-p(1));
if (A == 1) adr2 = ones(lf,1); end
if (A > 1) adr2 = [ ones(lf,1), adr ]; end
for (k = 2:n)
	% Applying each IFS mapping:
	f2 = [ f2, p(k) + phi(k)*(f-p(k)) ];
	if (A == 1) adr2 = [ adr2; k*ones(lf,1) ]; end
	if (A > 1) adr2 = [ adr2; k*ones(lf,1), adr ]; end
end

f = f2;
if (A >= 1) adr=adr2; else adr=0; end

end
