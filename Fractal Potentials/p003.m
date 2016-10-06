
% Plot fractal flow stream lines.

clf;
hold on;
clear all;

% Inputs:
optfra = input('Fractal [1..4] = ');
N = input('Level [quick=5, medium=9, accurate=12] = ');
visopt = input('Type of plot [LC=1, VF=2, SL=3] = '); % 1 = level curves; 2 = vector field; 3 = stream lines

% Secondary inputs:
offs = 0.2; % offset percentage
numdiv = 1000; % number of coordinate partitions
LC = 150; % number of level values for the level curves
stabopt = 1; % stabilize?
vis = 1; % 1 = rectangular, 2 = ellipsoidal, 3 = focal visualization

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

if (visopt == 1)
	
	% Plot the level curves:
	divval = 0:(2*pi/LC):2*pi;
	contourf(x,y,psival,divval)
	% contour3(x,y,psival,divval)
	% contour(x,y,psival,divval)
else
	% Calculate vector field:
	[ u1, u2 ] = gradient(psival, dx, dy);
	if (stabopt == 1)
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

% Set axes:
axis equal;
axis(SetAxes(pts, offs*100));
axis off;

% Plot fixed points:
plot(real(p), imag(p), 'wo')

% Plot fractal:
plot(real(pts), imag(pts), 'y.')

hold off;
