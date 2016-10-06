% Enclosing circle radius for points p from a fixed center c.
% Input:	p, c
% Output:	radius

function res = Rho(p, c)

res = max(abs(p-c));
