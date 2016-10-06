% Plot a bifractal with its convex hull.

clear all;

% Inputs:
N = 12; % fractal level, quick=9, accurate=12
mrsz = 3; % pt size, 3 or 10

% Set the fractal parameters:
flp = 1;
demopt = input('Demo parameters [1..4] = ');
switch (demopt)
case 1
	phi = [ (0.5^flp)*exp((2*flp*pi*0)*i), (sqrt(0.5)^flp)*exp((2*flp*pi/8)*i) ]
	p = [ 1+0.5*i, 0+1*i ]
case 2
	phi = [ (0.7^flp)*exp((-pi*flp/6)*i), (0.59^flp)*exp((pi*flp/4)*i) ]
	p = [ 0 1 ]
case 3
	phi = [ (0.7^flp)*exp((pi*flp/3)*i), (0.6^flp)*exp((pi*flp/5)*i) ]
	p = [ 0 1 ]
case 4
	phi = [ (0.65^flp)*exp((-2*flp*pi/6)*i), (0.65^flp)*exp((2*flp*pi/4)*i) ]
	p = [ 0 1 ]
end

% Plot the fractal:
clf;
hold on;
pts = IFS(N, phi, p, p);
plot(real(pts), imag(pts), 'k.', 'Markersize', mrsz)

% Plot main fixed points and their first iterates:
plot(real(p), imag(p), 'ro')
T1p2 = p(1)+phi(1)*(p(2)-p(1));
T2p1 = p(2)+phi(2)*(p(1)-p(2));
plot(real(T1p2), imag(T1p2), 'ro')
plot(real(T2p1), imag(T2p1), 'ro')

% Convex hull:
[chull, e] = FCHull2(phi, p);
% plot(real(chull), imag(chull), 'ro-')
plot(real(e), imag(e), 'bo')
t = -0.3:0.01:20;
e1t = p(1)+(phi(1).^t).*(e(1)-p(1));
e2t = p(2)+(phi(2).^t).*(e(2)-p(2));
plot(real(e1t), imag(e1t), 'c-')
plot(real(e2t), imag(e2t), 'c-')

% Set axes:
axis equal;
axis(SetAxes([pts e1t e2t], 5));
axis off;

hold off;
