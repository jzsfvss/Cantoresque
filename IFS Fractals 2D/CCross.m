% Calculates the "cross product" of two complex numbers.

function res = CCross(p1,p2)

res = real(p1)*imag(p2)-imag(p1)*real(p2);