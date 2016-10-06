% Plot the  t -> |u + e^ct v|_c  curve.

clear all;

% Input:
tminfac = 1; % >=1
tmaxfac = 2; % >=1
offs = 0.5; % offset
urad = 5; % in [1,10]
vrad = 1; % in [1,10]

% Initializing parameters:
xmax = 10;
ymax = 10;
uvi = input('u,v [1] random, [2] demo, [3] input = ');
switch (uvi)
case 1
	u = urad*(rand+i*rand);
	v = vrad*(rand+i*rand);
	phi = rand*exp(i*pi*(rand-1)); % Arg in (-pi, 0]
case 2
	u = 1.2755+2.5298*i;
	v = 0.6991+0.8909*i;
%	phi = 0.5*exp(-i*2*pi/3);
	phi = 0.5*exp(-i*pi/3);
otherwise
	u = input('u = ');
	v = input('v = ');
	phi = input('phi = ');
end

% Calculating the spiral value:
tau = (Arg(u)-Arg(v))/Arg(phi);
tdiv = pi/abs(Arg(phi));
tmin = tau-tminfac*tdiv;
tmax = tau+tmaxfac*tdiv;
t = (tmin-offs):0.001:(tmax+offs);
tlen = length(t);
st = Spi(u+(phi.^t)*v, phi, 0);
su = Spi(u,phi,0);
ymin = 0.9*min(st);
ymax = 1.1*max(st);
y = ymin:0.01:ymax;
ylen = length(y);

% Plotting:
clf;
hold on;
plot(t, st, 'b-');
plot(t, Spi(u,phi,0)*ones(1,tlen), 'r--');
for tck = tmin:tdiv:tmax
	plot(tck*ones(1,ylen), y, 'r--');
end
plot((tau-0.5)*ones(1,ylen), y, 'y--');
plot((tau+0.5)*ones(1,ylen), y, 'y--');
plot(round(tau)*ones(1,ylen), y, 'g--');
plot((round(tau)-1)*ones(1,ylen), y, 'g--');
plot((round(tau)+1)*ones(1,ylen), y, 'g--');
% axis equal;
axis([ tmin, tmax, ymin, su+(su-ymin) ]);
% set(gca, 'XTick', tmin:tdiv:tmax);
% set(gca, 'YTick', su);
axis off;
hold off;
