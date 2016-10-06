
% Plot a fractal.

clear all;

% Inputs:
[ phi, p ] = IFSdm();
N = 4; % fractal level, quick=9, accurate=12
mrsz = 1; % pt size, 3 or 10

% Set the fractal parameters:
dopts = 1; % plot fixed points? 
dotit = 0; % 1/0=y/n - title?
doptlab = 0; % 1/0=y/n - pt labels?
% phi = [ (0.5^flp)*exp((2*flp*pi*0)*i), (sqrt(0.5)^flp)*exp((2*flp*pi/8)*i) ]; % Fractal 1
% phi = [ 0.7*exp((-pi/6)*i), 0.59*exp((pi/4)*i) ]; % Fractal 2
% phi = [ 0.7*exp((pi/3)*i), 0.6*exp((pi/5)*i) ]; % Fractal 3
% phi = [ (0.65^flp)*exp((-2*flp*pi/6)*i), (0.65^flp)*exp((2*flp*pi/4)*i) ]; % Fractal 4
% p = [ 0 1 ];
% p = [ 1+0.5*i, 0+1*i ];
S = p;

% Plot the fractal:
pts = IFS(N, phi, p, S);
clf;
hold on;
plot(real(pts), imag(pts), 'k.', 'Markersize', mrsz)

% Set axes:
axis equal;
axis(SetAxes(pts, 5));
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
