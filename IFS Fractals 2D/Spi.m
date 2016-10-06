% Spiral absolute value.

function sv = Spi(z, phi, p)

b = Pitch(phi);

sv = abs(z-p)./exp(b*Arg(z-p));
