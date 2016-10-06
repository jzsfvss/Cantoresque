% Convex hull of bifractals (theoretically).

function [chull, e] = FCHull2(phi, p)

e = ExtPts2(phi, p, 50); % cycle of the princ extr pt
pts = e;

% Map it around:
for k = 1:2
	Nk = ceil(abs(2*pi/angle(phi(k))));
	Tne = e;
	for n = 1:Nk
		Tne = p(k)+phi(k)*(Tne-p(k));
		pts = [ pts, Tne ];
	end
end


chind = convhull(real(pts), imag(pts));
chull = pts(chind);
