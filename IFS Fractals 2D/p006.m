
% Plot a 2-map fractal, with the general and the circumcircle.

clear all;
clf;
hold on;

% Inputs:
N = 12; % fractal level, quick=9, accurate=12
mrsz = 3; % pt size, 3 or 10

% Set the fractal parameters:
flp = 1; % power of the flow
doptlab = 0; % 1/0=y/n - pt labels?

% Example 1:
% phi = [ (sqrt(0.5)^flp)*exp((2*flp*pi/8)*i), (0.5^flp)*exp((2*flp*pi*0)*i) ];
% p = [ 1+0.5*i, 0+1*i ];
% Example 2:
phi = [ (0.65^flp)*exp((-2*flp*pi/6)*i), (0.65^flp)*exp((2*flp*pi/4)*i) ];
p = [ 0 1 ];
% Example 3:
% phi = [ 0.7*exp((pi/3)*i), 0.6*exp((pi/5)*i) ];
% p = [ 0 1 ];

S = p;

% Plot the fractal:
pts = IFS(N, phi, p, S);
plot(real(pts), imag(pts), 'k.', 'Markersize', mrsz)

% Plot main fixed points:
offs = 0.03;
plot(real(p), imag(p), 'ro')
if (doptlab)
	text(real(p(1))+offs, imag(p(1))-offs, 'p1')
	text(real(p(2))+offs, imag(p(2))-offs, 'p2')
end
T1p2 = p(1)+phi(1)*(p(2)-p(1));
T2p1 = p(2)+phi(2)*(p(1)-p(2));
plot(real(T1p2), imag(T1p2), 'ro')
plot(real(T2p1), imag(T2p1), 'ro')
if (doptlab)
	text(real(T1p2)+offs, imag(T1p2)-offs, 'T1(p2)')
	text(real(T2p1)+offs, imag(T2p1)-offs, 'T2(p1)')
end

% Circumcircle:
cc = CircleCir(phi, p);
ccp = Circle(cc(1), cc(2));
plot(real(ccp), imag(ccp), 'r')
pts = [pts, ccp];

% General circle:
cg = CircleGen(phi, p);
cgp = Circle(cg(1), cg(2));
plot(real(cgp), imag(cgp), 'b')
pts = [pts, cgp];

% Set axes:
axis equal;
axis(SetAxes(pts, 5));
axis off;

% Generating the title:
title(['circumcircle radius = ', num2str(cc(2)), ';  general circle radius = ', num2str(cg(2)) ])

hold off;
