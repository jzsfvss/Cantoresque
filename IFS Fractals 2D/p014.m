% Convex hull for bifractals (numerically)

clear all;

% Input:
[ phi, p ] = IFSdm();
L = 20; % level
A = 20; % address length
mrsz = 1; % pt size 1-10
pemth = input('Principal extremal pt (1) numerically (2) theoretically (3) guess = ');

% Fractal points:
pts = IFS(L, phi, p, p(1));

% Extremal points:
[e, adre, pe1, pe2] = ExtPtsNM(phi, p, L, A, pemth);

% Convex hull:
chind = convhull(real(e), imag(e));
e = e(chind);

% Plotting the fractal:
clf;
hold on;
T1pts = p(1)+phi(1)*(pts-p(1));
T2pts = p(2)+phi(2)*(pts-p(2));
plot(real(T1pts), imag(T1pts), 'r.', 'Markersize', mrsz)
plot(real(T2pts), imag(T2pts), 'b.', 'Markersize', mrsz)
plot(real(p), imag(p), 'ro')
T1p = p(1)+phi(1)*(p-p(1));
T2p = p(2)+phi(2)*(p-p(2));
plot(real(T1p), imag(T1p), 'k.', 'Markersize', 3)
plot(real(T2p), imag(T2p), 'k.', 'Markersize', 3)

% Plotting the convex hull and its first iterates:
plot(real(e), imag(e), 'ro-')
T1e = p(1)+phi(1)*(e-p(1));
T2e = p(2)+phi(2)*(e-p(2));
plot(real(T1e), imag(T1e), 'k-')
plot(real(T2e), imag(T2e), 'k-')
plot(real(T1e), imag(T1e), 'm.', 'Markersize', 2)
plot(real(T2e), imag(T2e), 'm.', 'Markersize', 2)

% Set axes:
axis equal;
axis(SetAxes(pts, 5));
axis off;

hold off;
