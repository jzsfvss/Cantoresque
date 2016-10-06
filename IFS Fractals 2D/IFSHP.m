% Returns the IFS iterates of S to a certain level N (high precision).
% N = 		iteration level (natural number)
% phi = 	factors of contraction (row vector of complex numbers)
% p = 		fixed points (row vector of complex numbers)
% S = 		initial set (row vector of complex numbers)
% dsplev =	display level (0/1)

function f = IFSHP(N, phi, p, S, dsplev)

n=length(phi);

if (N == 0)

f = S;
w = S*0+1;

else

if (dsplev), disp([ num2str(N), ' more levels to go...' ]); end

% Previous level points:
f = IFSHP(N-1, phi, p, S, dsplev);

f2 = mp(p(1)+phi(1)*(f-p(1)));
for (k = 2:n)
	% Applying each IFS mapping:
	f2 = [ f2, mp(p(k)+phi(k)*(f-p(k))) ];
end

f = f2;

end
