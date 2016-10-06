% Set title for rational angled fractals.
% Input:	phi,p
% Output:	title string vector

function res = SetTitRat(phi, p)

sp = ' ';
res1 = '[lk] [thk] [pk] =';
res2 = mat2str(abs(phi));
res3 = rats(angle(phi)/(2*pi), 6);
res4 = mat2str(p);

res = [ res1 sp res2 sp '[' res3 ']' sp res4 ];
