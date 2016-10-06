% Extremal points theoretically.

clear all;

% Input:
[ phi, p, demopt ] = IFSdm();
L = 20; % level
A = 500; % address length

% Fractal points:
f = IFS(L, phi, p, p(1));

% Extremal points:
[ epr, ecyc, eall, adrpr, adrcyc, adrall ] = ExtPts2(phi, p, A);
size(adrall)

% Plotting:
clf;
hold on;
plot(real(f), imag(f), 'k.', 'Markersize', 3);
plot(real([eall, eall(1)]), imag([eall, eall(1)]), 'ro-');
plot(real(ecyc), imag(ecyc), 'mo');
plot(real(epr), imag(epr), 'm*');
plot([0 1], [0 0], 'co');
for k = 1:length(eall)
	tloc = (1+0.1/abs(eall(k)))*eall(k);
	%text(real(tloc), imag(tloc), num2str(k));
end
if (demopt == 5)
	n1 = abs(round(angle(phi(1))*360/pi));
	n2 = round(angle(phi(2))*360/pi);
	l1 = round(abs(phi(1))*100)/100;
	l2 = round(abs(phi(2))*100)/100;
	title([ num2str(l1) 'exp(-i' num2str(n1) 'pi/' num2str(360) '),  ' num2str(l2) 'exp(i' num2str(n2) 'pi/' num2str(360) ')' ]);
else
	title([ num2str(abs(phi(1))) 'exp(' num2str(angle(phi(1))) 'i), ' num2str(abs(phi(2))) 'exp(' num2str(angle(phi(2))) 'i)' ]);
end
disp([[1:length(eall)]', conj(eall')]);
AdrConvDisp(adrall);
axis equal;
axis(SetAxes([f, eall, tloc ], 5));
axis off;
hold off;

% Figure:
svfig = input('Save figure? 1/0 = ');
if (svfig)
	export_fig tmp -png -transparent -m3;
end
