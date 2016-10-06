% Neighbouring indices in a matrix.

function neivec = MatNeiInd(row,col,rows,cols)

neir = 0;
neic = 0;

for r = (row-1):(row+1)
	for c = (col-1):(col+1)
		if (r >= 1 && r <= rows && c >= 1 && c <= cols && (sum([r,c]==[row,col]) ~= 2))
			neir = [neir, r];
			neic = [neic, c];
		end
	end
end
neir = neir(2:end);
neic = neic(2:end);

neivec = sort(sub2ind([rows cols], neir, neic));
