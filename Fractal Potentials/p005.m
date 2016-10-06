
% Plot fractal flow stream lines.

clear all;

% Inputs:
optfra = input('Fractal [1..4] = ');
N = input('Level [quick=5, medium=9, accurate=12] = ');
rar = input('Streamline density [dense=10, medium=50, rare=100] = ');
pltgrd = input('Plot grid? [y/n]=[1/0] = '); % plot grid?

% Secondary inputs:
stabthr = 2; % threshold for stabilization [1.1-10]
offs = 0.2; % offset percentage
numdiv = 1000; % number of coordinate partitions
dostab = 1; % do stabilization?
colmin = [0 0 1]; % blue
colmax = [0 1 1]; % cyan
colnum = round(100/rar); % number of colors

% IFS:
switch (optfra)
case 1
	phi = [ 0.5*exp((2*pi*0)*i), sqrt(0.5)*exp((2*pi/8)*i) ]; % Fractal 1
	p = [ 1+0.5*i, 0+1*i ];
case 2
	phi = [ 0.7*exp((-pi/6)*i), 0.59*exp((pi/4)*i) ]; % Fractal 2
	p = [ 0, 1 ];
case 3
	phi = [ 0.7*exp((pi/3)*i), 0.6*exp((pi/5)*i) ]; % Fractal 3
	p = [ 0, 1 ];
otherwise
	phi = [ 0.65*exp((-2*pi/6)*i), 0.65*exp((2*pi/4)*i) ]; % Fractal 4
	p = [ 0, 1 ];
end
prb = [ 0.5, 0.5 ]; % probabilities

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
miny = miny-offs*diamy;
maxy = maxy+offs*diamy;
dx = (maxx-minx)/numdiv;
dy = (maxy-miny)/numdiv;

% Generate points:
x = minx:dx:maxx;
y = miny:dy:maxy;
lx = length(x);
ly = length(y);
Z = zeros(ly,lx);
psival = zeros(ly,lx);
for k = 1:lx
	Z(:,k) = x(k)+y*i;
end

% Calculate stream values:
disp('Calculating stream values...')
psival = PsiIFS(Z, phi, p, prb, N);

% Calculate vector field:
disp('Calculating gradient field...')
[ u1, u2 ] = gradient(psival, dx, dy);

% Stabilization:
if (dostab)

disp('Stabilization of gradient field...')
u = u1 + u2*i;
nu = u;
au = abs(u);
cnt = 0;
aveau = mean(mean(au));
minau = min(min(au));
maxau = max(max(au));
%disp(['Mean gradient norm = ', num2str(aveau)]);
%disp(['Min gradient norm = ', num2str(minau)]);
%disp(['Max gradient norm = ', num2str(maxau)]);
for k = 1:ly
	for j = 1:lx
		neicur = MatNeiInd(k,j,ly,lx);
		cn = 0;
		sn = 0;
		ln = length(neicur);
		mnn = min(au(neicur));
		mxn = max(au(neicur));
		if (mxn > stabthr*mnn)
			for n = 1:ln
				if (au(neicur(n)) < stabthr*mnn)
					sn = sn + u(neicur(n));
					cn = cn+1;
				end
			end
			nu(k,j) = sn/cn;
			cnt = cnt+1;
		end
	end
	if (mod(k,numdiv/10) == 0)
		disp(['Stabilized ', num2str(100*k/numdiv), '%']);
	end
	% disp(['Made ', num2str(cnt), ' corrections in row ', num2str(k), '/', num2str(ly), '.']);
end
u1 = -imag(nu);
u2 = real(nu);

else

nu1 = u1;
nu2 = u2;
u1 = -nu2;
u2 = nu1;

end

% Starting points for streamlines:
disp('Grid for streamlines...')
rarx = x(1:rar:lx);
lrx = length(rarx);
stx = rarx;
sty = miny+zeros(1,lrx);
for k = 1:(numdiv/rar)
	stx = [stx, rarx];
	sty = [sty, miny + k*rar*dy + zeros(1,lrx)];
end
lstx = length(stx);
lsty = length(sty);

clf;
hold on;

% Plot grid:
if (pltgrd)
	plot(stx,sty,'y.');
	axis equal;
	axis([minx maxx miny maxy]);
end

% Streamlines:
% streamline(x,y,u1,u2,stx,sty);
colvec = colmax-colmin;
for k=1:colnum
	disp(['Calculating streamlines - Stage ', num2str(k), '/', num2str(colnum), '...'])
	h = streamline(x, y, u1, u2, stx(k:colnum:lstx), sty(k:colnum:lstx));
	set(h, 'Color', colmin+(k/colnum)*colvec);
	view(2);
end

% Set axes:
axis equal;
axis([minx maxx miny maxy]);
axis off;

% Plot fractal:
plot(real(p), imag(p), 'ro')
plot(real(pts), imag(pts), 'r.')
disp('Streamlines and fractal plotted...')

sv = input('Save figure? [no/jpg/png/both]=[0/1/2/3] ');
if (sv == 1 || sv == 3)
	savefig('tmp','jpeg');
	disp('JPEG saved.')
end
if (sv == 2 || sv == 3)
	export_fig tmp -png -transparent -r300;
	disp('PNG saved.')
end

hold off;
