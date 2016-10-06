% Plot an evolving Sierpinski triangle

clear all;

% Input:
L = 4; % max iteration level
sca = 3; % scaling factor
pltfull = 1; % plot the attractor?
%filcol = 'r';
filcol = 0.7*[ 1, 1, 1 ]; % grey = color of the fill
phi = [ 0.5, 0.5, 0.5 ];
p = [ 0, 1, (1+sqrt(3)*i)/2 ];
S = 0;

% Generate and plot;
clf;
hold on;
fill(sca*real(p), sca*imag(p), filcol, 'LineWidth', 0.2);
for l = 1:L
	%display(l)
	pts = IFS(l, phi, p, S);
	N = length(pts);
	mrk = 1.125*l;
	for n = 1:N
		pn = sca*(mrk + pts(n) + p/(2^l));
		fill(real(pn), imag(pn), filcol, 'LineWidth', 0.2);
	end
end
if (pltfull)
	pts = IFS(12, phi, p, S);
	mrk = 1.125*L + 1.5;
	plot(sca*(mrk + real(pts)), sca*imag(pts), 'k.', 'Markersize', 1);
	%plot(sca*(mrk + [real(p), 0]), sca*[imag(p), 0], 'k-', 'LineWidth', 0.2);
	plot(sca*(mrk - [ 0.4, 0.25, 0.1 ]), sca*[ 1, 1, 1 ]*sqrt(3)/4, '.', 'Markersize', 9, 'Color', filcol/2);
	xmax = sca*(1.125*L + 2.5);
else
	xmax = sca*(1.125*L + 1);
end

% Set axes:
axis equal;
axis(SetAxes([ sca*p, xmax ], 20));
axis off;
hold off;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
axis equal;

% Save;
sv = input('Save? 0/1 = ');
%{
if (sv)
	set(gcf, 'Color', 'w');
	export_fig tmp -jpg -m6 -painters;
	display('Plot saved to tmp.jpg.');
end
%}
if (sv)
	addpath('..\..\Open Source\export_fig\');
	export_fig p00066e_Sierpinski -png -transparent -m3;
	%addpath('..\..\Open Source\plot2svg\');
	%plot2svg('p00066e_Sierpinski.svg');
end