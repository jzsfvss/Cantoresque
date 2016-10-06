% Extremal points numerically (no brainer).
% L = 	iteration level (natural number)
% A = 	max length of address
% phi = factors of contraction (row vector of complex numbers)
% p = 	fixed points (row vector of complex numbers)
% pemth = primary extrema (1) numerically or (2) theoretically

function [e, adre, pe1, pe2, kpe1, kpe2] = ExtPtsNM(phi, p, L, A, pemth)

tolt = 0.001; % translational tolerance
tola = 0.001; % angular tolerance

% Generate the fractal and the addresses:
[ f, adrf ] = IFSA(phi, p, L, A);

% Find the convex hull extrema:
ind = convhull(real(f), imag(f))';
li = length(ind);
k = 1;
while (k <= (li-2))
	if (abs(f(ind(k))-f(ind(k+1))) < tolt)
		ind = [ind(1:k), ind((k+2):li)];
		li = li-1;
	end
	k = k+1;
end
k = 2;
while (k <= (li-1))
	cond = (abs(angle((f(ind(k-1))-f(ind(k)))/(f(ind(k))-f(ind(k+1))))) < tola);
	if (cond)
		ind = [ind(1:(k-1)), ind((k+1):li)];
		li = li-1;
	end
	k = k+1;
end
if (abs(f(ind(li-1))-f(ind(li))) < tolt || abs(f(ind(1))-f(ind(li))) < tolt)
	ind = ind(1:(li-1));
	li = li-1;
end

% Assign the extrema:
e = f(ind);
adre = adrf(ind,:);

% Primary extrema:
if (pemth)
%	s1 = mod(log(Spi(e, phi(1), p(1))), 2*pi);
	s1 = Spi(e, phi(1), p(1));
	if (angle(phi(1)) < 0)
		ind1 = find(s1 == max(s1));
	else
		ind1 = find(s1 == min(s1));
	end
	pe1 = e(ind1);
	s2 = mod(log(Spi(e, phi(2), p(2))), 2*pi);
	ind2 = find(s2 == max(s2));
	pe2 = e(ind2);
else
	[pe adr] = ExtPts2(phi, p, A);
	pe1 = pe(1);
	pe2 = pe(2);
	ind1 = find(abs(e-pe1) < 0.01);
	ind2 = find(abs(e-pe2) < 0.01);
	ind2 = ind2(1);
end
k = 1;
while (adre(ind1,k) == 1) k=k+1; end
pe1 = AdrB(p(1), [adre(ind1, k:A); ones(1,A-k+1)], phi, p);
mpe1 = min(abs(e-pe1));
pe1 = e(find(abs(e-pe1) < mpe1+0.0001));
if (length(pe1) > 1) pe1=pe1(1); end
k = 1;
while (adre(ind2,k) == 2) k=k+1; end
pe2 = AdrB(p(1), [adre(ind2, k:A); ones(1,A-k+1)], phi, p);
mpe2 = min(abs(e-pe2));
pe2 = e(find(abs(e-pe2) < mpe2+0.0001));
if (length(pe2) > 1) pe2=pe2(1); end

k1 = 0;
k2 = 0;
kpe1 = find(abs(e-pe1) < 0.0001);
kpe1 = kpe1(1);
kpe2 = find(abs(e-pe2) < 0.0001);
kpe2 = kpe2(1);
