% Returns the IFS iterates of S to a certain level N along with the weights.
% N = 	iteration level (natural number)
% phi = factors of contraction (row vector of complex numbers)
% p = 	fixed points (row vector of complex numbers)
% pr = 	probabilities (row vector of [0,1] numbers)
% S = 	initial set (row vector of complex numbers)

function [f, w] = IFSW(N, phi, p, pr, S)

n=length(phi);

if (N == 0)

f = S;
w = S*0+1;

else

% Previous level points:
[f, w] = IFSW(N-1, phi, p, pr, S);

f2 = p(1) + phi(1)*(f-p(1));
w2 = pr(1)*w;
for (k = 2:n)
	% Applying each IFS mapping:
	f2 = [ f2, p(k) + phi(k)*(f-p(k)) ];
	w2 = [ w2, pr(k)*w ];
end

f = f2;
w = w2;

end
