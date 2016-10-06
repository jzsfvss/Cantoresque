% Sets the axes to a certain set of complex points in the plane.
% Creates a percentage offset around the points (can be 0).
% Input: 	pts, offperc
% Output: 	[xmin xmax ymin ymax].

function res = SetAxes(pts, offperc)

xmin = min(real(pts));
xmax = max(real(pts));
ymin = min(imag(pts));
ymax = max(imag(pts));

xoffset = (offperc/100)*(xmax-xmin);
yoffset = (offperc/100)*(ymax-ymin);

xmin = xmin-xoffset;
xmax = xmax+xoffset;
ymin = ymin-yoffset;
ymax = ymax+yoffset;

res = [xmin xmax ymin ymax];