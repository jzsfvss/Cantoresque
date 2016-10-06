addpath('..\savefig\');
addpath('..\export_fig\');
addpath('..\Multiprecision Computing Toolbox\');

if input('Change digits to 50? 0/1 = ')
	mp.Digits(70); % 34 = quadruple precision computation, optimized efficiency
end
if input('Show longest number format? 0/1 = ')
	format longG % force MatLab to show all digits of numbers
end

% mp.Test() % tests each program
