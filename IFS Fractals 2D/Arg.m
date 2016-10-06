function a = Arg(c)

a = angle(c);
badind = find(a == -pi);
a(badind) = pi;