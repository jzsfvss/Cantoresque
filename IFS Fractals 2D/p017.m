% Plot a bifractal with the circumcircle iterated to some level.

clear all;

% Inputs:
[ phi, p ] = IFSdm();
N = 20; % fractal iteration level, quick=9, accurate=12
L = input('Circle iteration level = ');
highres = 1;
plotlabels = 0;
if (highres)
	mrsz = 0.25; % pt size, 1 or 3 or 10
	ccmrsz = 2;
	clw = 0.2; % Circle linewidth, default = 0.5
else
	mrsz = 1;
	ccmrsz = 3;
	clw = 0.3;
end
if (plotlabels)
	docc = 1; % plot circumcircle?
	dopts = 1; % plot fixed points? 
	dotit = 1; % 1/0=y/n - title?
	doptlab = 1; % 1/0=y/n - pt labels?
else
	docc = 0; % plot circumcircle?
	dopts = 0; % plot fixed points? 
	dotit = 0; % 1/0=y/n - title?
	doptlab = 0; % 1/0=y/n - pt labels?
end

% Generating the fractal:
S = p;
fpts = IFS(N, phi, p, S);

% Generating the circumcircle iterates:
cc = CircleCir(phi, p);
c = cc(1);
r = cc(2);
[cit, rit] = IFSW(L, phi, p, abs(phi), c);
rit = rit*r;

% Plotting:
clf;
hold on;
plot(real(fpts), imag(fpts), 'k.', 'Markersize', mrsz)
ccpts = 0;
if (docc)
	ccpts = Circle(c, r, 0.01);
	plot(real(ccpts), imag(ccpts), 'r', 'LineWidth', clw)
	plot(real(c), imag(c), 'r.', 'Markersize', ccmrsz)
end
cirpts = 0;
for (k = 1:length(cit))
	cirk = Circle(cit(k),rit(k), 0.01);
	plot(real(cirk), imag(cirk), 'r', 'LineWidth', clw)
	plot(real(cit), imag(cit), 'r.', 'Markersize', ccmrsz)
	cirpts = [cirpts, cirk];
end

% Set axes:
axis equal;
axis(SetAxes([fpts, cirpts, ccpts], 10));
axis off;

% Plot main fixed points:
offs = 0.03;
if (dopts)
	plot(real(p), imag(p), 'ro')
end
if (doptlab)
	text(real(p(1))+offs, imag(p(1))-offs, 'p1')
	text(real(p(2))+offs, imag(p(2))-offs, 'p2')
end
if (dopts)
	T1p2 = p(1)+phi(1)*(p(2)-p(1));
	T2p1 = p(2)+phi(2)*(p(1)-p(2));
	plot(real(T1p2), imag(T1p2), 'yo')
	plot(real(T2p1), imag(T2p1), 'yo')
end
if (doptlab)
	text(real(T1p2)+offs, imag(T1p2)-offs, 'T1(p2)')
	text(real(T2p1)+offs, imag(T2p1)-offs, 'T2(p1)')
end

% Generating the title:
if (dotit) title(SetTitRat(phi, p)); end

hold off;
