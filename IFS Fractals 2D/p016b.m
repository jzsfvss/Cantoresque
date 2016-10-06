% Convex hull for bifractals (simple rounding).

clear all;

% Base input:
L = 20; % fractal iteration level
N = 180; % collected address length for optimization (>20)
mrsz = 1; % pt size 1-10

% Input:
[ phi, p, demopt, angrat ] = IFSdm();

% Fractal points:
pts = IFS(L, phi, p, p(1));

% Extremal points:
[ ecyc, eall, adrpr ] = ExtPtsRnd(phi, p, N);
per = length(ecyc);
disp([ 'Principal extremum = ', num2str(ecyc(1)) ])
disp([ 'Principal address = ', AdrConvStr(AdrConv(adrpr)) ])

% Testing if ecyc(1) is extremal:
%chind = convhull(real(pts), imag(pts));
%enm = pts(chind); % numerical extrema
%isext = ~inpolygon(real(ecyc(1)), imag(ecyc(1)), real(enm), imag(enm));
%if (isext)
%	disp('It is indeed extremal.')
%else
%	disp('It is not extremal.')
%end

% Plotting the fractal:
clf;
hold on;
plot(real(pts), imag(pts), 'k.', 'Markersize', mrsz);
T1p = p(1)+phi(1)*(p(2)-p(1));
T2p = p(2)+phi(2)*(p(1)-p(2));
plot(real(p), imag(p), 'ro')
plot(real(T1p), imag(T1p), 'ro')
plot(real(T2p), imag(T2p), 'ro')
plot(real(p), imag(p), 'r.')

% Plotting the convex hull:
plot(real(eall), imag(eall), 'r-')
plot(real(ecyc(2:per)), imag(ecyc(2:per)), 'ro')
plot(real(ecyc(1)), imag(ecyc(1)), 'r*')

% Set axes:
axis equal;
axis(SetAxes(pts, 5));
axis off;

% Title:
title([ '\lambda_1=', num2str(abs(phi(1))), ', \lambda_2=', num2str(abs(phi(2))), ', \theta_1=', num2str(angrat(1,1)), '\pi/', num2str(angrat(1,2)), ', \theta_2=', num2str(angrat(2,1)), '\pi/', num2str(angrat(2,2)), ' (rounding)' ]);

hold off;
