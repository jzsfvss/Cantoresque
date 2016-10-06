% Returns the IFS iterates of a set S to a certain level N.
% N = 	iteration level (natural number)
% phi = factors of contraction (row vector of complex numbers)
% p = 	fixed points (row vector of complex numbers)
% S = 	initial set (row vector of complex numbers)

function res = IFS(N, phi, p, S)

n=length(phi);

if (N == 1)

res = S;

else

% Previous level points:
res = IFS(N-1, phi, p, S);

res2 = p(1) + phi(1)*(res - p(1));
for (k = 2:n)
	% Applying each IFS mapping:
	res2 = [ res2, p(k) + phi(k)*(res - p(k)) ];
end

res = [ res, res2 ];

end