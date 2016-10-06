% Plot a 3-map fractal with the circumcircle.

clear all;

% Inputs:
N = input('Level [12] = '); % fractal level, quick=9, accurate=12
mrsz = 1; % pt size, 1, 3, 10
demopt = input('Demo parameters [1..7] = ');

switch (demopt)
case 1
	phi = [ 0.5, 0.5, 0.5 ]
	p = [ 0, 1, 0.3+0.7*i ]
case 2
	phi = [ 0.5*exp((2*pi/24)*i), 0.5*exp((2*pi/48)*i), 0.5*exp((2*pi/32)*i) ]
	p = [ 0, 1, 0.3+0.7*i ]
case 3
	phi = [ 0.5*exp((2*pi/12)*i), 0.5*exp((2*pi/24)*i), 0.5*exp((2*pi/16)*i) ]
	p = [ 0, 1, 0.3+0.7*i ]
case 4
	phi = [ sqrt(0.5)*exp((2*pi/8)*i), 0.5*exp((2*pi*0)*i), 0.7*exp((2*pi*2)*i) ]
	p = [ 1+0.5*i, 0+1*i, 0 ]
case 5
	phi = [ 0.65*exp((-2*pi/6)*i), 0.65*exp((2*pi/4)*i), 0.65*exp((2*pi/5)*i) ]
	p = [ 0, 1, 0.3+0.7*i ]
case 6
	phi = [ 0.7*exp((pi/3)*i), 0.6*exp((pi/5)*i), 0.5*exp((pi/7)*i) ]
	p = [ 0, 1, 0.3+0.5*i ]
case 7
	phi = [ 0.65*exp((-2*pi/6)*i), 0.65*exp((2*pi/4)*i), 0 ]
	p = [ 0, 1, 0.3+0.7*i ]
end

%pr = abs(phi).^DimSim(phi); % OLD

% Plot the fractal:
clf;
hold on;
S = p(1);
%[pts, w] = IFSW(N, phi, p, pr, S); % OLD
pts = IFS(N, phi, p, S);
plot(real(pts), imag(pts), 'k.', 'Markersize', mrsz)

%v1=length(pts)
%v2=length(w)

% Plot first maps of fixed points:
Hp = IFS(1, phi, p, p);
plot(real(Hp), imag(Hp), 'ro')
% Plot main fixed points:
plot(real(p), imag(p), 'ro')

% Circumcircle:
cc = CircleCir(phi, p)
ccp = Circle(cc(1), cc(2), 0.01);
for (k = 1:3)
	cck = [p(k)+phi(k)*(cc(1)-p(k)), abs(phi(k))*cc(2)];
	cckp = Circle(cck(1), cck(2), 0.01);
	plot(real(cckp), imag(cckp), 'r--')
end
plot(real(ccp), imag(ccp), 'r')
pts = [pts, ccp];

% Set axes:
axis equal;
axis(SetAxes(pts, 5));
axis off;
set(gcf, 'Color', 'w'); % set background to white

hold off;
