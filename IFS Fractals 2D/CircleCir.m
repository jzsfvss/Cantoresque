% Calculates the center and radius of the circumcircle for 2-map or 3-map IFS.
% Input: 	phi, p
% Output:	[c, r]

function res = CircleCir(phi, p)

if (length(p) == 2)

ph1 = phi(1);
ph2 = phi(2);
p1 = p(1);
p2 = p(2);
l1 = abs(ph1);
l2 = abs(ph2);
mu = (l2-l1)/(1-0.5*(l1+l2));

c = ((2-mu)*(1-ph1)*p1 + (2+mu)*(1-ph2)*p2)/(4 - (2-mu)*ph1 - (2+mu)*ph2);
r = 2*abs((1-ph1)*(1-ph2)*(p2-p1))/((1-0.5*(l1+l2))*abs(4-(2-mu)*ph1-(2+mu)*ph2));

res = [c, r];

end

if (length(p) == 3)

alp2 = ((1-abs(phi))./abs(1-phi)).^2;
ap2 = abs(p).^2;

A = (alp2(3)-alp2(2))*p(1) + (alp2(1)-alp2(3))*p(2) + (alp2(2)-alp2(1))*p(3);
B = 2*CCross(p(2)-p(1),p(2)-p(3));
C = (ap2(2)-ap2(3))*p(1) + (ap2(3)-ap2(1))*p(2) + (ap2(1)-ap2(2))*p(3);
c0 = C/(B*i);
r0 = abs(c0-p(1));
c1 = A/(B*i);
D = (1/B)*(alp2(1)*CCross(p(2),p(3)) + alp2(2)*CCross(p(3),p(1)) + alp2(3)*CCross(p(1),p(2))) + CDot(c0,c1);

if (abs(c1) < 0.00001)
	c = c0;
	r = r0;
else
	discr = (D^2)-abs(c1^2)*(r0^2);
	if (discr < 0)
		c = 0;
		r = 0;
	else
		r2 = (-D - sqrt(discr))/(abs(c1^2));
		r = sqrt(r2);
		c = c0+c1*(r2);
	end
end

res = [c, r];

end
