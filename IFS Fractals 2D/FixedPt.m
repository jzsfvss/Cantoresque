% Calculate fixed point with a given address.
% adr is in two-row vector format.
% Bimap IFS only.

function fp = FixedPt(adr, phi, p)

fp = AdrB(0, adr, phi, p)/(1 - Phi(adr, phi)); % Based on the Slope Lemma
