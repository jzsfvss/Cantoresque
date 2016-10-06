
% Plot a 3-map fractal with the circumcircle.

clear all;

% Inputs:
L = input('Level [9/12] = '); % fractal level, quick=9, accurate=12
mrsz = 3; % pt size, 3 or 10
demopt = input('Demo parameters [1..7] = ');
vdir = input('Angle of light rays [0..180] = ');
v = exp(vdir*pi*i/180)/i; % normal to the light rays
prec = input('Precision factor [1..10] = ');

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
	phi = [ 0.65*exp((-2*pi/6)*i), 0.65*exp((2*pi/4)*i) ]
	p = [ 0, 1 ]
end

pr = abs(phi).^DimSim(phi); % natural similarity probability weights

% Plot the fractal:
clf;
hold on;
[f, w] = IFSW(L, phi, p, pr, p(1));
plot(real(f), imag(f), 'k.', 'Markersize', mrsz)

% Calculate density function:
z0 = min(CDot(f,v))*v;
zN = max(CDot(f,v))*v;
eps0 = (min(abs(phi))^(L/prec))*abs(p(2)-p(1));
N = round(abs(zN-z0)/eps0);
eps = abs(zN-z0)/N; % spacing between light rays
h1 = min(CDot(f,v*i));
h2 = max(CDot(f,v*i));
endpts = [z0 zN]+h1*v*i;
plot(real(endpts), imag(endpts), 'b-');
z = z0+(0:N)*eps*v;
ypr = 0*z;
for j = 1:length(f)
	k = 1+ceil(CDot(f(j)-z0,v)/eps);
	if (k > length(ypr)) k = length(ypr); end
	ypr(k) = ypr(k)+w(j);
end
y = z + ((ypr/max(ypr))*(h2-h1)+h1)*v*i;
plot(real(y), imag(y), 'r-');

% Set axes:
axis equal;
axis(SetAxes([f, y, endpts], 5));
axis off;

hold off;
