% Extremal points numerically (no brainer).

clear all;
diary
delete('tmp.txt')
diary('tmp.txt')

% Input:
[ phi, p ] = IFSdm();
L = 12; % level
A = 20; % address length
pemth = input('Primary extrema (1) numerically (2) theoretically (3) guess = ');

% Fractal points:
f = IFS(L, phi, p, 1.6896-0.6491i);
% f = IFS(L, phi, p, p(1));

% Extremal points:
if (pemth < 3)
	[e, adre, pe1, pe2, kpe1, kpe2] = ExtPtsNM(phi, p, L, A, pemth);
else
	[e, adre, pe1, pe2] = ExtPtsNM(phi, p, L, A, 1);
	kpe1 = input('Primary extremum wrt T1 = ');
	pe1 = e(kpe1);
	kpe2 = input('Primary extremum wrt T2 = ');
	pe2 = e(kpe2);
end

% Plotting:
clf;
hold on;
disp(' ');
plot(real(f), imag(f), 'k.', 'Markersize', 1);
plot(real([e, e(1)]), imag([e, e(1)]), 'ro-');
plot(real([pe1 pe2]), imag([pe1 pe2]), 'bo');
t = 0:0.01:20;
e1t = p(1)+(phi(1).^t)*(pe1-p(1));
e2t = p(2)+(phi(2).^t)*(pe2-p(2));
plot(real(e1t), imag(e1t), 'm-', 'LineWidth', 0.25)
plot(real(e2t), imag(e2t), 'c-', 'LineWidth', 0.25)
for k = 1:length(e)
	tloc = (1+0.1/abs(e(k)))*e(k);
	text(real(tloc), imag(tloc), num2str(k));
end
disp([[1:length(e)]', conj(e')]);
AdrConvDisp(adre);
disp(' ')
disp(['Pr extr wrt T1:  e(', num2str(kpe1), ') = ', num2str(pe1)]);
AdrConvDisp(adre(kpe1,:));
disp(['Pr extr wrt T2:  e(', num2str(kpe2), ') = ', num2str(pe2)]);
AdrConvDisp(adre(kpe2,:));
axis equal;
%axis(SetAxes([f, tloc], 5));
axis(SetAxes([f, tloc, e1t, e2t], 5));
axis off;
hold off;

diary off