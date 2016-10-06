% Load demo parameters.
% angrat is a 2x2 matrix with first column the numerator and second column the denominator of the angles (*pi).
% If the angle denominator is unknown then 360 is picked.

function [ phi, p, demopt, angrat, highprec ] = IFSdm()

global phi
global p

demopt = input('Input [0] / Random [1] / Demo [2..8] = ');
flp = input('Flow degree [<=1] = ');
highprec = input('High precision calculations & plotting? 0/1 = ');
angrat = zeros(2,2);

switch (demopt)
case 0
	if (highprec)
		disp('Put apostrophes around the number...');
	end
	l1 = input('Lambda 1 = ');
	l2 = input('Lambda 2 = ');
	a1 = input('Angle 1 = ');
	a2 = input('Angle 2 = ');
	if (highprec)
		l1 = mp(l1);
		l2 = mp(l2);
		a1 = mp(a1);
		a2 = mp(a2);
		phi = [ mp((l1^flp)*exp(flp*a1*i)), mp((l2^flp)*exp(flp*a2*i)) ];
	else
		phi = [ (l1^flp)*exp(flp*a1*i), (l2^flp)*exp(flp*a2*i) ];
	end
	p = [ 0 1 ];
	angrat = [ double(round(180*a1/pi)), 180; double(round(180*a2/pi)), 180 ];
case 1
	disp('Random C-type bifractal:') % P <= Q, a1n <= a2n
	a2n = round(rand(1)*180);
	a1n = round(rand(1)*a2n);
	a1 = -pi*a1n/180;
	a2 = pi*a2n/180;
	l1 = round(rand(1)*10000)/10000;
	l2 = round(rand(1)*10000)/10000;
	if (highprec)
		phi = [ mp(l1*exp(a1*i)), mp(l2*exp(a2*i)) ]
	else
		phi = [ l1*exp(a1*i), l2*exp(a2*i) ]
	end
	p = [ 0 1 ]
	angrat = [ -a1n, 180; a2n, 180 ];
case 2
	disp('Fractal 1:')
	phi = [ (0.5^flp)*exp((2*flp*pi*0)*i), (sqrt(0.5)^flp)*exp((2*flp*pi/8)*i) ]
	p = [ 0+1*i, 1+0.5*i ]
	angrat = [ 0, 4; 1*flp, 4 ];
case 3
	disp('Fractal 2:')
	phi = [ (0.7^flp)*exp((-pi*flp/6)*i), (0.59^flp)*exp((pi*flp/4)*i) ]
	p = [ 0 1 ]
	angrat = [ -1*flp, 6; 1*flp, 4 ];
case 4
	disp('Fractal 3:')
	phi = [ (0.7^flp)*exp((5/30)*2*pi*i*flp), (0.6^flp)*exp((3/30)*2*pi*i*flp) ]
	p = [ 0 1 ]
	angrat = [ 5*flp, 15; 3*flp, 15 ];
case 5
	disp('Fractal 4:')
	phi = [ (0.65^flp)*exp((-2*flp*pi/6)*i), (0.65^flp)*exp((2*flp*pi/4)*i) ]
	p = [ 0 1 ]
	angrat = [ -2*flp, 6; 3*flp, 6 ];
case 6
	disp('Lévy C Curve:')
	phi = [ (sqrt(0.5)^flp)*exp((-1/8)*flp*2*pi*i), (sqrt(0.5)^flp)*exp((1/8)*flp*2*pi*i) ]
	p = [ 0 1 ]
	angrat = [ -1*flp, 4; 1*flp, 4 ];
case 7
	disp('Heighway Dragon:')
	phi = [ (sqrt(0.5)^flp)*exp((1/8)*flp*2*pi*i), (sqrt(0.5)^flp)*exp((3/8)*flp*2*pi*i) ]
	p = [ 0 1 ]
	angrat = [ 1*flp, 4; 3*flp, 4 ];
case 8
	disp('Twindragon / Davis-Knuth Dragon:')
	phi = [ (sqrt(0.5)^flp)*exp((-1/8)*flp*2*pi*i), (sqrt(0.5)^flp)*exp((3/8)*flp*2*pi*i) ]
	p = [ 0 1 ]
	angrat = [ -1*flp, 4; 3*flp, 4 ];
otherwise
	phi = [ (0.9^flp)*exp((-3*flp)*i), (0.5^flp)*exp((0.4*flp)*i) ]
	p = [ 0 1 ]
end
