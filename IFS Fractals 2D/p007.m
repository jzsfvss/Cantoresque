
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

% Circumcircle:
cc = CircleCir(phi, p);
ccp = Circle(cc(1), cc(2));
plot(real(ccp), imag(ccp), 'r')
pts = [pts, ccp];

% Set axes:
axis equal;
axis(SetAxes(pts, 5));
axis off;

% Plot a streamline of T1:
t = (-100):0.1:10;
ang = 0:0.1:2*pi;
lang = length(ang);
Na = length(ang);
for j = (lang-7)
	pts = p(1) + (phi(1).^t)*exp(ang(j)*i)*(p(2)-p(1));
	plot(real(pts), imag(pts), 'm');
end

hold off;
