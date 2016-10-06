% Creates circle points for plotting.
% Input:	c (complex), r, dth (~0.01)
% Output:	set of points

function res = Circle(c, r, dth)

res = c + r*exp([0:dth:2*pi+dth]*i);
