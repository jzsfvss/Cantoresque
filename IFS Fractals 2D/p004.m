
% Plot a vortex with iterates.

clear all;
clf;
hold on;

% Set vortex parameters:
phi = sqrt(0.5)*exp((2*pi/12)*i);
p = 0;

% Plot streamlines of T(z)=p+phi*(z-p):
t = (-100):0.1:10;
ang = 0:0.1:2*pi;
Na = length(ang);
for j = 1:Na
	pts = p + (phi.^t)*exp(ang(j)*i)*(1-p);
	plot(real(pts), imag(pts), 'b');
end

% Set axes:
axis equal;
axis([ -0.3, 1, -0.2, 0.5 ]);
axis off;

% Plot points:
plot(real(p), imag(p), 'ro')
for j = 1:4:Na
	z=exp(ang(j)*i);
	for k=1:30
		plot(real(z), imag(z), 'ro')
		z = p+phi*(z-p);
	end
end

hold off;
