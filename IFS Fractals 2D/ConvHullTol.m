% Convex hull with tolerance.
% e =	extrema
% ind = indices among the original points

function [ e, ind ] = ConvHullTol(f, tol)

tolt = tol; % translational tolerance
tola = tol; % angular tolerance

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

e = f(ind);
