% Convex hull for bifractals (linear optimization).

clear all;

% Base input:
L = 20;
% N = 1000; % collected address length for optimization (>20, recommended: 180) % use J below!
mrsz = 1; % pt size 1-10
dospi = 0; % draw spirals
dsplev = 0; % display generation levels
plaux = 1; % plot fixed pts and iterates

% Input:
[ phi, p, demopt, angrat, highprec ] = IFSdm();
cycper = (-angrat(1,1)+angrat(2,1))/gcd(-angrat(1,1),angrat(2,1)); % period predicted by theory
J = -angrat(1,1)/gcd(-angrat(1,1),angrat(2,1)); % collected period predicted by theory % use below!
L0 = input('Number of iterations = '); % fractal iteration level (typical: 20)
dofrac = input('Generate the fractal? 0/1 = ');
if (highprec && dofrac), dsplev = 1; end
dohull = input('Determine the convex hull? 0/1 = ');
if (dohull)
	docomp = input('Do Lin vs. Rnd comparison? 0/1 = ');
	if (dofrac)
		genst = input('Generate fractal from p_1 or e_*? 1/2 = ');
	end
else
	docomp = 0;
	genst = 1;
end

% Extremal points:
if (dohull)
	disp('Calculating the extremal points via linear optimization...')
	if (highprec)
		[ ecyc, eall, adrpr, useful ] = ExtPtsLinHP(phi, p, N);
	else
		[ ecyc, eall, adrpr, useful ] = ExtPtsLin(phi, p, J);
	end
	if (docomp)
		disp('Calculating the extremal points via simple rounding...')
		[ ecyc2, eall2, adrpr2 ] = ExtPtsRnd(phi, p, N);
	end
	per = length(ecyc);
	disp(' ')
	disp([ 'Principal extremal point = ', num2str(ecyc(1)) ])
	disp([ 'Principal address = ', AdrConvStr2(AdrConv(adrpr)) ])
end

% Testing if ecyc(1) beats ecyc2(1):
if (docomp)
	disp(' ')
	if (useful)
		uf='yes';
		disp('Linear optimization result differs from simple rounded result.')
	else
		uf='no';
		disp('Linear optimization result does not differ from simple rounded result.')
	end
	disp('Testing Lin vs. Rnd procedure...')
	beatsLinRnd = ~inpolygon(real(ecyc(1)), imag(ecyc(1)), real(eall2), imag(eall2));
	beatsRndLin = ~inpolygon(real(ecyc2(1)), imag(ecyc2(1)), real(eall), imag(eall));
	if (beatsLinRnd || beatsRndLin)
		if (beatsLinRnd)
			disp('Linear optimization beats the simple rounded procedure.')
			bts = 'Lin>Rnd';
		end
		if (beatsRndLin)
			disp('The simple rounded procedure beats linear optimization.')
			bts = 'Rnd>Lin';
		end
	else
		disp('Both results are on the convex hull.')
		bts = 'Lin=Rnd';
	end
else
	bts = 'NoComp';
end

% Generating the fractal points:
if (dofrac)
	if (dsplev)
		disp(['Generating the fractal up to level ', num2str(L), '...'])
	else
		disp('Generating the fractal...')
	end
	if (genst)
		if (highprec)
			pts = IFSHP(L, phi, p, p(1), dsplev);
		else
			pts = IFS(L, phi, p, p(1));
			% pts = IFS2(phi, p, 10, 3, 6000);
		end
	else
		if (highprec)
			pts = IFSHP(L, phi, p, ecyc(1), dsplev);
		else
			pts = IFS(L, phi, p, ecyc(1));
		end
		% numpts = length(pts)
	end
	if (L0 < 20)
		pts0 = IFS(L0, phi, p, p(1));
	end
end

% Plotting the fractal:
if (dofrac)
	disp('Plotting the fractal...')
	clf;
	hold on;
	plot(real(pts), imag(pts), 'k.', 'Markersize', mrsz);
	if (L0 < L)
		plot(real(pts0), imag(pts0), 'b.', 'Markersize', 3);
	end
	if (plaux)
		T1p = p(1)+phi(1)*(p(2)-p(1));
		T2p = p(2)+phi(2)*(p(1)-p(2));
		plot(real(p), imag(p), 'ro')
		plot(real(T1p), imag(T1p), 'mo')
		plot(real(T2p), imag(T2p), 'mo')
		plot(real(p), imag(p), 'r.')
		plot(real(T1p), imag(T1p), 'm.')
		plot(real(T2p), imag(T2p), 'm.')
	end
end

% Plotting the convex hull:
if (dohull)
	plot(real(eall), imag(eall), 'r-')
	plot(real(ecyc(1:per)), imag(ecyc(1:per)), 'ro')
	if (docomp), plot(real(ecyc2(1)), imag(ecyc2(1)), 'g.'), end
end

% Plot optimal line:
if (dohull)
	plot(real(ecyc(1)), imag(ecyc(1)), 'b.')
	ln0 = -1:0.01:1;
	v0 = (1-phi(2))*log(phi(1))/abs((1-phi(2))*log(phi(1)));
	ln = ecyc(1)+v0*ln0;
	plot(real(ln), imag(ln), 'b-')
end

% Plot spirals:
if (dospi)
	t = (-10):0.1:100;
	spi = p(1) + (phi(1).^t)*(ecyc(1)-p(1));
	plot(real(spi), imag(spi), 'g-');
	spi = p(2) + (phi(2).^t)*(ecyc(1)-p(2));
	plot(real(spi), imag(spi), 'g-');
end

% Set axes:
axis equal;
if (dohull)
	axis(SetAxes([ pts, eall ], 5));
else
	axis(SetAxes(pts, 5));
end
axis off;

% Labels:
mytit = [ '\lambda_1=', num2str(double(round(abs(phi(1))*10000)/10000)), ', \lambda_2=', num2str(double(round(abs(phi(2))*10000)/10000)), ', \theta_1=', num2str(angrat(1,1)), '\pi/', num2str(angrat(1,2)), ', \theta_2=', num2str(angrat(2,1)), '\pi/', num2str(angrat(2,2)), ', L=', num2str(L) ];
if (dohull), mytit = [ mytit, ', period=', num2str(cycper), ', ', bts ]; end
if (highprec), mytit = [ mytit, ', high prec' ]; end
title(mytit);

hold off;
