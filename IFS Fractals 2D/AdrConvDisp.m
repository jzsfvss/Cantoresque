% Displays a matrix of addresses as strings.

function res = AdrConvDisp(mat)

sz = size(mat);
A = sz(2);
N = sz(1);

if (N == 1)
	disp([ 'Address:  ', AdrConvStr(AdrConv(mat(1,1:A))) ]);
else
	for (k = 1:N)
		if (k < 10) sp=' '; else sp=''; end
		disp([ sp num2str(k) ':  ' AdrConvStr(AdrConv(mat(k,1:A))) ]);
	end
end

res = 1;
