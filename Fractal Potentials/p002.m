
% Plot fractal stream lines.

clear all;
clf;

% Input:
numdiv = 10*100;
N = 5; % fractal level, quick=9, accurate=12
LC = 150; % number of level values for the level curves
[ phi, p ] = IFSdm();
% phi = [ 0.7*exp((pi/3)*i), 0.6*exp((pi/5)*i) ]; % factors
% p = [ 0, 1 ]; % fixed pt
prb = [ 0.5, 0.5 ]; % probabilities
offs = 0.2; % offset percentage
opt = 1; % 1 = level curves; 2 = vector field; 3 = stream lines
stopt = 1; % stabilize?
vis = 3; % 1 = rectangular, 2 = ellipsoidal visualization

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
dx = 2*0.2*(maxx-minx)/numdiv;
dy = 2*0.2*(maxy-miny)/numdiv;
thr = 1+dx/10; % threshold for stabilization

% Generate points and psi values:
x = minx:dx:maxx;
y = miny:dy:maxy;
lx = length(x);
ly = length(y);
Z = zeros(ly,lx);
psival = zeros(ly,lx);
for k = 1:lx
	Z(:,k) = x(k)+y*i;
end
Z2 = reshape(Z',lx*ly,1); % makes a column vector of the pts

% Calculate stream values:
psival = PsiIFS(Z, phi, p, prb, N);

hold on;
if (opt == 1)
	% Plot the level curves:
	divval = 0:(2*pi/LC):2*pi;
	contourf(x,y,psival,divval)
	% contour3(x,y,psival,divval)
	% contour(x,y,psival,divval)
else
	% Calculate vector field:
	[ u1, u2 ] = gradient(psival, dx, dy);
	if (stopt == 1)
		nu1 = zeros(ly,lx);
		nu2 = zeros(ly,lx);
		for k = 1:ly
			for j = 1:lx
				nu1(k,j) = u1(k,j);
				nu2(k,j) = u2(k,j);
				if (k > 1 && abs((u1(k,j)+u2(k,j)*i)/(u1(k-1,j)+u2(k-1,j)*i)) > thr)
						nu1(k,j) = u1(k-1,j);
						nu2(k,j) = u2(k-1,j);
				end
				if (k < ly && abs((u1(k,j)+u2(k,j)*i)/(u1(k+1,j)+u2(k+1,j)*i)) > thr)
						nu1(k,j) = u1(k+1,j);
						nu2(k,j) = u2(k+1,j);
				end						
				if (j > 1 && abs((u1(k,j)+u2(k,j)*i)/(u1(k,j-1)+u2(k,j-1)*i)) > thr)
						nu1(k,j) = u1(k,j-1);
						nu2(k,j) = u2(k,j-1);
				end
				if (j < lx && abs((u1(k,j)+u2(k,j)*i)/(u1(k,j+1)+u2(k,j+1)*i)) > thr)
						nu1(k,j) = u1(k,j+1);
						nu2(k,j) = u2(k,j+1);
				end
				if ((k > 1) && (j > 1) && (abs((u1(k,j)+u2(k,j)*i)/(u1(k-1,j)+u2(k-1,j)*i)) > thr) && (abs((u1(k,j)+u2(k,j)*i)/(u1(k,j-1)+u2(k,j-1)*i)) > thr))
						nu1(k,j) = (u1(k-1,j) + u1(k,j-1))/2;
						nu2(k,j) = (u2(k-1,j) + u2(k,j-1))/2;
				end
			end
		end
		u1 = -nu2;
		u2 = nu1;
	end
	if (opt == 2)
		quiver(x,y,u1,u2);
	else
		switch (vis)
		case 1
			stx = [ minx+zeros(1,ly), x, maxx+zeros(1,ly), x ];
			sty = [ y, miny+zeros(1,lx), y, maxy+zeros(1,lx) ];
			lstx = length(stx);
			lsty = length(sty);
			stx = stx(1:(numdiv/100):lstx);
			sty = sty(1:(numdiv/100):lsty);
		case 2
			th1 = 0:0.2:2*pi;
			th2 = 0:0.1:2*pi;
			th3 = 0:0.05:2*pi;
			th4 = 0:0.05:2*pi;
			stx = [ 0.5+1*cos(th1), 0.5+0.75*cos(th2), 0.5+0.5*cos(th3), 0.5+0.25*cos(th4) ];
			sty = [ 1*sin(th1), 0.75*sin(th2), 0.5*sin(th3), 0.25*sin(th4) ];
		case 3
			th = 0:0.1:2*pi;
			pts = IFS(N, phi, p, p);
			cir = 0.25*exp(th*i);
			st = 0.5+4*cir;
			for (k = 1:length(pts))
				st = [ st, pts(k)+cir ];
			end
			stx = real(st);
			sty = imag(st);
		end
		streamline(x,y,u1,u2,stx,sty);
	end
end

% Plot fractal:
plot(real(pts), imag(pts), 'y.')
axis equal;
axis([minx maxx miny maxy]);
axis off;

hold off;
