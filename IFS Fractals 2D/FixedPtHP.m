% Calculate fixed point with a given address (high precision).
% adr is in two-row vector format.
% Bimap IFS only.

function fp = FixedPtHP(adr, phi, p)

fp = mp(AdrBHP(0, adr, phi, p)/(1 - PhiHP(adr, phi))); % Based on the Slope Lemma
