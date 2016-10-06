addpath('C:\Users\Angi es Jozsi\Documents\MATLAB\savefig\');
addpath('C:\Users\Angi es Jozsi\Documents\MATLAB\export_fig\');
addpath('C:\Users\Angi es Jozsi\Documents\MATLAB\Multiprecision Computing Toolbox\');

if input('Change digits to 50? 0/1 = ')
	mp.Digits(70); % 34 = quadruple precision computation, optimized efficiency
end
if input('Show longest number format? 0/1 = ')
	format longG % force MatLab to show all digits of numbers
end

% mp.Test() % tests each program
