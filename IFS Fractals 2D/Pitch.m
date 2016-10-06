% Pitch of a log spiral.

function b = Pitch(phi)

b = log(abs(phi))./Arg(phi);
