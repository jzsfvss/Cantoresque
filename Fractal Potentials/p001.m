
% Plot vortex stream lines.

clear all;
clf;
hold on;

% Plot a vortex:
phi = sqrt(0.5)*exp((2*pi/12)*i);
p = 0; % fixed pt
rp = 1; % ref pt

% Generate points and psi values:
x = -10:0.1:10;
y = -10:0.1:10;
lx = length(x);
ly = length(y);
pts = zeros(ly,lx);
z = zeros(ly,lx);
for k = 1:lx
	pts(:,k) = x(k)+y*i;
	z(:,k) = PsiVor(pts(:,k), phi, p, rp);
end
pts = reshape(pts',lx*ly,1); % makes a column vector of the pts

% Plot the level curves:
contourf(x,y,z)

% Set axes:
axis equal;
axis(SetAxes(pts, 5));
axis off;

hold off;
