% Principal extremal point for rational C-type bifractals (linear optimization, high precision).
% ----------------------------------------------------------------------------------------------
% phi =		factors
% p =		fixed points
% N =		max collected adr length (>20 can be large)
% Need:		th1 = -2*pi*P/M, th2 = 2*pi*Q/M, 0<P<=Q
% ----------------------------------------------------------------------------------------------
% ecyc(1) =	principal extremal point
% adrpr =	address of the pr ext pt 	(row address)
% ecyc =	cyclical extremal pts 		(incl epr)
% eall =	all ext pts 			(sorted angularly)
% ----------------------------------------------------------------------------------------------

function [ ecyc, eall, adrpr, useful ] = ExtPtsLinHP(phi, p, N)

% Initializing:
l1 = mp(abs(phi(1)));
l2 = mp(abs(phi(2)));
th1 = mp(Arg(phi(1)));
th2 = mp(Arg(phi(2)));
alp = mp(atan(log(l1)/th1)); % alpha
r = mp(-th2/th1); % angle ratio
adr = zeros(2,2*N); % collected address mtx

% Calculating epr and adrpr=adr(epr) via linear optimization:
s = 0;
useful = 0;
for j = 1:N
	s0 = s;
	s1 = floor(r*j);
	s2 = ceil(r*j);
	tau1 = mp((l1^s1)*cos(th1*s1 + th2*j + alp));
	tau2 = mp((l1^s2)*cos(th1*s2 + th2*j + alp));
	if (tau1 >= tau2), s = s1; else s = s2; end;
	n = s-s0;
	adr(1,2*j-1) = 2; % T2
	adr(2,2*j-1) = 1; % T2 power
	adr(1,2*j) = 1; % T1
	adr(2,2*j) = n; % T1 power
	useful = max(useful, s~=round(r*j) && abs(r*j-round(r*j))~=0.5); % Is doing the new linear optimization useful?
end
adrpr = AdrConv(adr); % convert to row adr of princ ext pt
L = length(adrpr); % len of adr of princ ext pt
[ per, ap ] = SeqPeriod(adrpr); % the period and cyclical part of adrpr
adrpr = ap;
epr = FixedPtHP(AdrConv(ap), phi, p);

% Calculating the cyclical addresses and extrema:
adrcyc = AdrCyc(ap); % row-by-row cyclic permutations of ap
ecyc = zeros(1,per);
ecyc(1) = epr;
for k = 2:per
	ecyc(k) = FixedPtHP(AdrConv(adrcyc(k,:)), phi, p);
end

% Calculating all extrema:
eall = ecyc;
K = 2*ceil(pi/abs(angle(phi(1))))+1; % needed iteration number for 180 deg
% Iterating to the left:
for k = 1:K
	eall = [ eall, mp(p(1)+(phi(1)^k)*(ecyc-p(1))) ];
end
K = 2*ceil(pi/abs(angle(phi(2))))+1;  % needed iteration number for 180 deg
for k = 1:K
	eall = [ eall, mp(p(2)+(phi(2)^k)*(ecyc-p(2))) ];
end
% Keeping the relevant points:
chind = convhull(double(real(eall)), double(imag(eall)));
eall = eall(chind);
