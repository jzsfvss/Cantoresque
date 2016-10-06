
% Plot vortex stream lines based on the velocity field.

clear all;

colmin = [0 0 1]; % blue
colmax = [0 1 1]; % cyan

[ phi, p ] = IFSdm();

xmin = -10;
xmax = 10;
xdiv = 0.1;
xdiv2 = 0.5;
xvec = xmin:xdiv:xmax;
xvec2 = xmin:xdiv2:xmax;
xnum = length(xvec);
ymin = -10;
ymax = 10;
ydiv = 0.1;
ydiv2 = 0.5;
yvec = ymin:ydiv:ymax;
yvec2 = ymin:ydiv2:ymax;
ynum = length(yvec);
[ x, y ] = meshgrid(xvec, yvec);
z = x+y*i;
zb = [ xmin+yvec2*i, xvec2+ymin*i, xmax+yvec2*i, xvec2+ymax*i ];
zbl = length(zb);

v = VecVor(z, phi(1), p(1));
v = min(xdiv,ydiv)*v/max(max(abs(v)));

hold on;

for k=1:zbl
	h = streamline(real(z), imag(z), real(v), imag(v), real(zb(k)), imag(zb(k)));
	set(h, 'Color', colmin+rand(1)*(colmax-colmin));
	view(2);
end

t = 0:0.1:100;
s0 = xmin+yvec(round(ynum/2))*i;
s = p(1)+(phi(1).^t)*(s0-p(1));
plot(real(s), imag(s), 'r-');

% quiver(x, y, real(v), imag(v));

% Set axes:
axis equal;
axis([xmin xmax ymin ymax]);
axis off;

hold off;
