% Returns the randomized points of a bimap IFS fractal.
% L =	max total iteration level (~100)
% N = 	collected iteration level (~30)
% Nf = 	number of points to be generated (~4000)
% phi = factors of contraction (row vector of two complex numbers)
% p = 	fixed points (row vector of two complex numbers)

function f = IFS2(phi, p, L, N, Nf)

% Generate power matrix:
pmt = zeros(Nf,N+1);
%mxn = 2*(L-N)/N;
for r = 1:Nf
	for c = 2:(N+1)
		pmt(r,c) = pmt(r,c-1) + round(-0.5+rand*(1+mxn));
	end
end
pmt = phi(1).^pmt;

% Generate fractal points:
phi1p = phi(1).^round(-0.5+rand(Nf,1)*(1+mxn));
phi2p = phi(2).^(0:1:N)';
f = phi1p.*(pmt*phi2p);
f = (1-phi(2))*f;
f = p(1) + f*(p(2)-p(1));
