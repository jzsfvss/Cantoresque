
% Plot fractal flow stream lines.

clear all;

% Inputs:
optfra = input('Fractal [1..4] = ');
N = input('Level [quick=5, medium=9, accurate=12] = ');

% Secondary inputs:
offs = 0.2; % offset percentage
numdiv = 1000; % number of coordinate partitions
blnum = 1; % number of blocks per row / column
blwid = numdiv/blnum;
dostab = 1; % Do stabilization?

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
stabthr = 1+dx/10; % threshold for stabilization

% Generate points:
allx = minx:dx:maxx;
ally = miny:dy:maxy;
lallx = length(allx);
lally = length(ally);

for blrow = 1:blnum
for blcol = 1:blnum

clf;
hold on;

x = allx((blwid*(blcol-1)+1):(blwid*blcol));
y = ally((lally-blwid*blrow+1):(lally-blwid*(blrow-1)));
lx = length(x);
ly = length(y);
Z = zeros(ly,lx);
psival = zeros(ly,lx);
for k = 1:lx
	Z(:,k) = x(k)+y*i;
end

% Calculate stream values:
psival = PsiIFS(Z, phi, p, prb, N);

% Calculate vector field:
[ u1, u2 ] = gradient(psival, dx, dy);

% Stabilization:
if (dostab)

nu1 = zeros(ly,lx);
nu2 = zeros(ly,lx);
for k = 1:ly
	for j = 1:lx
		nu1(k,j) = u1(k,j);
		nu2(k,j) = u2(k,j);
		if (k > 1 && abs((u1(k,j)+u2(k,j)*i)/(u1(k-1,j)+u2(k-1,j)*i)) > stabthr)
			nu1(k,j) = u1(k-1,j);
			nu2(k,j) = u2(k-1,j);
		end
		if (k < ly && abs((u1(k,j)+u2(k,j)*i)/(u1(k+1,j)+u2(k+1,j)*i)) > stabthr)
			nu1(k,j) = u1(k+1,j);
			nu2(k,j) = u2(k+1,j);
		end						
		if (j > 1 && abs((u1(k,j)+u2(k,j)*i)/(u1(k,j-1)+u2(k,j-1)*i)) > stabthr)
			nu1(k,j) = u1(k,j-1);
			nu2(k,j) = u2(k,j-1);
		end
		if (j < lx && abs((u1(k,j)+u2(k,j)*i)/(u1(k,j+1)+u2(k,j+1)*i)) > stabthr)
			nu1(k,j) = u1(k,j+1);
			nu2(k,j) = u2(k,j+1);
		end
		if ((k > 1) && (j > 1) && (abs((u1(k,j)+u2(k,j)*i)/(u1(k-1,j)+u2(k-1,j)*i)) > stabthr) && (abs((u1(k,j)+u2(k,j)*i)/(u1(k,j-1)+u2(k,j-1)*i)) > stabthr))
			nu1(k,j) = (u1(k-1,j) + u1(k,j-1))/2;
			nu2(k,j) = (u2(k-1,j) + u2(k,j-1))/2;
		end
	end
end
u1 = -nu2;
u2 = nu1;

else

nu1 = u1;
nu2 = u2;
u1 = -nu2;
u2 = nu1;

end

%divval = 0:(2*pi/50):2*pi;
%contourf(x,y,psival,divval)
%pause(1);

stx = [ min(x)+zeros(1,ly), x, max(x)+zeros(1,ly), x ];
sty = [ y, min(y)+zeros(1,lx), y, max(y)+zeros(1,lx) ];
lstx = length(stx);
lsty = length(sty);
stx = stx(1:(numdiv/(100*1)):lstx);
sty = sty(1:(numdiv/(100*1)):lsty);

% streamline(x,y,u1,u2,stx,sty);
plot([min(x) max(x)], [min(y) max(y)], 'b.');
h=streamline(x,y,u1,u2,stx,sty);
set(h, 'Color', 'red');
view(2);

% Set axes:
axis equal;
axis([min(x) max(x) min(y) max(y)]);
axis off;

savefig('tmp','jpeg');
% print -djpeg 'tmp';
[myit, mp] = imread('tmp.jpeg');
if (blrow == 1 && blcol == 1)
	myits = size(myit);
	di = myits(1)-17;
	dj = myits(2)-19;
%	myi = zeros(blnum*di, blnum*dj);
end
if (blcol == 1)
	myi = zeros(di, blnum*dj);
end

myit = myit(16:(15+di), 4:(3+dj));
% myi((di*(blrow-1)+1):(di*blrow), (dj*(blcol-1)+1):(dj*blcol)) = myit;
myi(1:di, (dj*(blcol-1)+1):(dj*blcol)) = myit;

disp(['Block [', int2str(blrow), ',', int2str(blcol), '/', int2str(blnum), ',', int2str(blnum), '] rendered.']);

if (blcol == blnum)
	imwrite(myi, ['./Render/tmp', int2str(blrow), '.jpeg']);
	disp(['Row [', int2str(blrow), '/', int2str(blnum), '] rendered.']);
end
% imwrite(myi, ['./Render/tmp', int2str(blrow), int2str(blcol), '.jpeg']);

end
end

% imwrite(myi,'tmp.jpeg');

% Plot fractal:
%plot(real(p), imag(p), 'wo')
%plot(real(pts), imag(pts), 'y.')

hold off;
