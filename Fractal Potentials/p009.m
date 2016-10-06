% Plot fractal flow stream lines. Version 7. Based on the measure transform.

clear all;

% Inputs:
[ phi, p ] = IFSdm();
% phiconj = EdConj(phi);
ppr = p(1);
N = input('Level [quick=5, medium=9, accurate=12] = ');

% Secondary inputs:
prb = [ 0.5, 0.5 ]; % probabilities
offs = 0.3; % offset percentage
numdiv = 1000; % number of coordinate partitions
mrsz = 3; % pt size, 3 or 10
psidiv = 200; % number of level curve divisions
strwid = 0.1; % LineWidth of streamlines, default=0.5

% Fractal with weights:
disp('Generating the fractal...')
[ pts, w ] = IFSW(N, phi, p, prb, ppr);
minx = min(real(pts));
maxx = max(real(pts));
diamx = maxx-minx;
if (diamx == 0) diamx = 1; end
miny = min(imag(pts));
maxy = max(imag(pts));
diamy = maxy-miny;
if (diamy == 0) diamy = 1; end
minx = minx-offs*diamx;
maxx = maxx+offs*diamx;
if (miny == maxy)
	miny = miny-2*offs*diamx;
	maxy = maxy+2*offs*diamx;
else
	miny = miny-offs*diamy;
	maxy = maxy+offs*diamy;
end
dx = 2*(maxx-minx)/numdiv;
dy = 2*(maxy-miny)/numdiv;
xvec = minx:dx:maxx;
yvec = miny:dy:maxy;
lx = length(xvec);
ly = length(yvec);

% Stream values:
psival = zeros(ly,lx);
for j = 1:ly
	for k = 1:lx
		z = xvec(k)+yvec(j)*i;		
		psival(j,k) = w*mod(angle(z-pts),2*pi)';
	end
end
psival = psival/(2*pi);

clf;
hold on;

% Plot the level curves:
psimin = min(min(psival));
psimax = max(max(psival));
psid = (exp(psimax)-exp(psimin))/psidiv;
divexpval = exp(psimin):psid:exp(psimax);
goodind = find(log(divexpval) > -100);
divval = log(divexpval(goodind));
[c, h] = contourf(xvec, yvec, psival, divval);
set(h, 'Fill', 'off', 'LineWidth', strwid); 
view(2)

% Plot fractal:
% plot(real(p), imag(p), 'ro')
plot(real(pts), imag(pts), 'r.', 'Markersize', mrsz)
axis equal;
axis([minx maxx miny maxy]);
axis off;

% Rendering:
rend = input('Render? [0/1] ');
if (rend)
	export_fig tmp -png -transparent -m7;
end

hold off;
