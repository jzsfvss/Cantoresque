% Returns the IFS iterates of S to a certain level N.
% N = 		iteration level (natural number)
% phi = 	factors of contraction (row vector of complex numbers)
% p = 		fixed points (row vector of complex numbers)
% S = 		initial set (row vector of complex numbers)

function f = IFS(N, phi, p, S)

n=length(phi);

if (N == 0)

f = S;

else

% Previous level points:
f = IFS(N-1, phi, p, S);

f2 = p(1) + phi(1)*(f-p(1));
for (k = 2:n)
	% Applying each IFS mapping:
	f2 = [ f2, p(k) + phi(k)*(f-p(k)) ];
end

f = f2;

end
