% Plot fractal flow stream lines. Version 5.

clear all;

% Inputs:
[ phi, p ] = IFSdm();
N = input('Level [quick=5, medium=9, accurate=12] = ');
rarx = input('Streamline density [dense=1, rare=10] = ');
rary = rarx;

% Secondary inputs:
prb = [ 0.5, 0.5 ]; % probabilities
offs = 0.3; % offset percentage
strvrts = 17500; % number of curve vertices
numdiv = 1000; % number of coordinate partitions
numdiv2 = 1000/rarx; % number of streamline partitions
colmin = [0 0 1]; % blue
colmax = [0 1 1]; % cyan
colvec = colmax-colmin;
colnum = round(100/rarx); % number of colors
stgs = 10;
mrsz = 1; % pt size, 3 or 10

% Generate the fractal:
disp('Generating the fractal...')
pts = IFS(N, phi, p, p);
minx = min(real(pts));
maxx = max(real(pts));
diamx = maxx-minx;
miny = min(imag(pts));
maxy = max(imag(pts));
diamy = maxy-miny;
minx = minx-offs*diamx;
maxx = maxx+offs*diamx;
if (miny == maxy)
	miny = miny-2*offs*diamx;
	maxy = maxy+2*offs*diamx;
else
	miny = miny-offs*diamy;
	maxy = maxy+offs*diamy;
end
if (N <= 6)
	strwid = 0.5; % LineWidth of streamlines
else
	strwid = 0.2;
end
if (N <= 7)
	dx = 2*strwid*(maxx-minx)/numdiv;
	dy = 2*strwid*(maxy-miny)/numdiv;
else
	numdivx = round((maxx-minx)/((min(abs(phi))^N)*abs(p(2)-p(1))/5));
	numdivy = round((maxy-miny)/((min(abs(phi))^N)*abs(p(2)-p(1))/5));
	dx = (maxx-minx)/numdivx;
	dy = (maxy-miny)/numdivy;
	rarx = 0.5*(maxx-minx)/(1000*dx);
	rary = 0.5*(maxy-miny)/(1000*dy);
end
xvec = minx:dx:maxx;
xvec2 = minx:(rarx*dx):maxx;
yvec = miny:dy:maxy;
yvec2 = miny:(rary*dy):maxy;

% Generate grid:
[ x, y ] = meshgrid(xvec, yvec);
Z = x+y*i;

% Starting points for streamlines:
% Zb = [ (minx+2*rarx*dx)+yvec2*i, xvec2+(miny+2*rary*dy)*i, (maxx-2*rarx*dx)+yvec2*i, xvec2+(maxy-2*rary*dy)*i ];
Zb = [ minx+yvec2*i, xvec2+miny*i, maxx+yvec2*i, xvec2+maxy*i ];
Zbl = length(Zb);

% Calculate velocity field:
disp('Calculating the velocity field...')
vec = VecIFS(Z, phi, p, prb, N);
% vec = vec/max(max(abs(vec)));
u1 = real(vec);
u2 = imag(vec);

clf;
hold on;

% Streamlines:
tst = 1;
for j=0:(stgs-1)
	tic;
	indj = find(mod(1:Zbl, stgs) == j);
	for k=indj
		h = streamline(x, y, u1, u2, real(Zb(k)), imag(Zb(k)), [ 0.1, strvrts ]);
		set(h, 'Color', colmin+rand(1)*(colmax-colmin), 'LineWidth', strwid); % Default LineWidth=0.5
		view(2);
	end
	if (tst)
		disp(['Calculating streamlines for ~' num2str(ceil(toc*(stgs-1)/60)) ' mins...']);
		tst = 0;
	end
end

% Set axes:
axis equal;
axis([minx maxx miny maxy]);
axis off;

% Plot fractal:
% plot(real(p), imag(p), 'ro')
plot(real(pts), imag(pts), 'k.', 'Markersize', mrsz)
disp('Streamlines and fractal plotted...')
%clear all;
clearvars;

% Saving the figure:
sv = input('Save figure? [0/1] ');
if (sv)
	rsn = input('Size? [S/M/L]=[1/2/3] ');
	switch (rsn)
		case 3
			disp('Estimated time: ~18 mins...');
			export_fig tmp -png -transparent -m7;
		case 2
			disp('Estimated time: ~4 mins...');
			export_fig tmp -png -transparent -m4;
		otherwise
			disp('Estimated time: ~2 mins...');
			export_fig tmp -png -transparent -m2;
	end		
	disp('Image saved.')
end

hold off;
