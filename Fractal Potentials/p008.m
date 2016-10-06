% Plot fractal flow stream lines. Version 6. High resolution render in parts.

clear all;

% Inputs:
[ phi, p ] = IFSdm();
N = input('Level [quick=5, medium=9, accurate=12] = ');
typ = input('Type of the initial flow [superpos=1, char=2] = ');
harcon = input('Plot harmonic conjugate instead? [0/1] = ');

% Secondary inputs:
prb = [ 0.5, 0.5 ]; % probabilities
offs = 0.1; % offset percentage
% strvrts = 25000;
% strvrts = round(25000*(1+0.5*harcon));
strvrts = 1.5*round(25000*(1+0.5*harcon)); % streamline vertices
numdiv = 1000; % number of coordinate partitions
colmin = [0 0 1]; % blue
colmax = [0 1 1]; % cyan
colvec = colmax-colmin;
stgs = 7; % number of stages, default=7
stgst = 1; % stage start, default=1
mrsz = 1; % pt size, 3 or 10
strwid = 0.2; % LineWidth of streamlines, default=0.5

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
dx = 2*strwid*(maxx-minx)/numdiv;
dy = 2*strwid*(maxy-miny)/numdiv;
xvec = minx:dx:maxx;
yvec = miny:dy:maxy;

% Generate grid:
[ x, y ] = meshgrid(xvec, yvec);
Z = x+y*i;

% Starting points for streamlines:
% Zb = [ (minx+2*dx)+yvec*i, xvec+(miny+2*dy)*i, (maxx-2*dx)+yvec*i, xvec+(maxy-2*dy)*i ];
Zb = [ minx+yvec*i, xvec+miny*i, maxx+yvec*i, xvec+maxy*i ];
if (harcon)
	charifs = CharIFS(phi, prb);
	srcs = find(sign(real(i*sign(imag(charifs))*log(phi))) > 0);
	numsrcs = length(srcs);
	cir = (0.5*(dx+dy)/min(abs(phi(srcs)))^N)*exp(2*pi*(1:16)*i/16);
	srcscir = p(srcs(1))+cir;
	for k = 2:numsrcs
		srcscir = [ srcscir, p(srcs(k))+cir ];
	end
	ifssrcs = IFS(N, phi, p, srcscir);
	Zb = [ Zb, ifssrcs ];
	% plot(real(ifssrcs), imag(ifssrcs), 'r.');
	% pause
end
Zbl = length(Zb);

% Calculate velocity field:
lu = input('Load the velocity field in u.mat? [0/1] = ');
if (lu)
	load u.mat;
	disp('Velocity field loaded.')
else
	disp('Calculating the velocity field...')
	vec = VecIFS(Z, phi, p, prb, N, typ);
	% vec = vec/max(max(abs(vec)));
	if (harcon) vec = i*sign(imag(charifs))*vec; end
	u1 = real(vec);
	u2 = imag(vec);
	su = input('Save velocity field? [0/1] = ');
	if (su)
		save('u.mat', 'u1', 'u2')
		disp('Velocity field saved.')
		rq = input('Finished? [0/1] = ');
		if (rq) return; end
	end
end

% Compare velocity field u.mat to uoth.mat:
cu = input('Compare u.mat to uoth.mat? [0/1] = ');
if (cu)
	u11 = u1;
	u12 = u2;
	load uoth.mat;
	u21 = u1;
	u22 = u2;
	diffu = min(max(1, sqrt(((u11-u21).^2)+((u12-u22).^2))), 3);
    diffuave = mean(mean(diffu));
%	diffu = max(0, log(diffu));
	maxdiff = max(max(diffu));
	divval = log(1:(maxdiff-1)/10:maxdiff);
	[c, h] = contourf(xvec,yvec,log(diffu),log(divval));
	set(h, 'LineColor', 'none');
	axis equal;
	axis([minx maxx miny maxy]);
	axis off;
%	export_fig tmpcon -png -transparent -m5;
	rq = input('Finished? [0/1] = ');
	if (rq) return; end
	u1 = u11;
	u2 = u12;
end

% Streamlines:
disp('Generating streamlines and saving partial images...')
tst = 1;
for j=(stgst-1):(stgs-1)

clf;
hold on;

tic;
indj = find(mod(1:Zbl, stgs) == j);
for k=indj
	h = streamline(x, y, u1, u2, real(Zb(k)), imag(Zb(k)), [ 0.1, strvrts ]);
	set(h, 'Color', colmin+rand(1)*(colmax-colmin), 'LineWidth', strwid); 
	view(2);
end


% Plot fractal:
% plot(real(p), imag(p), 'ro')
plot(real(pts), imag(pts), 'k.', 'Markersize', mrsz)
axis equal;
axis([minx maxx miny maxy]);
axis off;

% Rendering:
%export_fig tmp -png -transparent -m3; % low res
export_fig tmp -png -transparent -m7;
if (j < 9)
	movefile('tmp.png', ['pic_0', num2str(j+1), '.png']);
else
	movefile('tmp.png', ['pic_', num2str(j+1), '.png']);
end

% Time estimation:
if (tst)
	disp(['Generating streamlines for ~' num2str(ceil(toc*(stgs-stgst)/60)) ' mins...']);
	tst = 0;
end

disp(['Generation stage ', num2str(j+1), '/', num2str(stgs), ' complete.'])
hold off;

% if (harcon) break; end

end % ending j-loop of streamlines

disp('Generation complete.')
