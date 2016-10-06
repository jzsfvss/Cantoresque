% Principal extremal point for rational C-type bifractals (with simple rounding).
% ------------------------------------------------------------------------------------------
% phi =		factors
% p =		fixed points
% N =		max collected adr length (>20 can be large)
% Need:		th1 = -2*pi*P/M, th2 = 2*pi*Q/M, 0<P<=Q
% ------------------------------------------------------------------------------------------
% ecyc(1) =	principal extremal point
% adrpr =	address of the pr ext pt 	(row address)
% ecyc =	cyclical extremal pts 		(incl epr)
% eall =	all ext pts 			(sorted angularly)
% ------------------------------------------------------------------------------------------

function [ ecyc, eall, adrpr ] = ExtPtsRnd(phi, p, N)

% Initializing:
l1 = abs(phi(1));
l2 = abs(phi(2));
th1 = Arg(phi(1));
th2 = Arg(phi(2));
alp = atan(log(l1)/th1); % alpha
r = -th2/th1; % angle ratio
adr = zeros(2,2*N); % collected address mtx

% Calculating epr and adrpr=adr(epr) via linear optimization:
s = 0;
for j = 1:N
	s0 = s;
	s = round(r*j);
	n = s-s0;
	adr(1,2*j-1) = 2; % T2
	adr(2,2*j-1) = 1; % T2 power
	adr(1,2*j) = 1; % T1
	adr(2,2*j) = n; % T1 power
end
adrpr = AdrConv(adr); % convert to row adr of princ ext pt
L = length(adrpr); % len of adr of princ ext pt
[ per, ap ] = SeqPeriod(adrpr); % the period and cyclical part of adrpr
adrpr = ap;
epr = FixedPt(AdrConv(ap), phi, p);

% Calculating the cyclical addresses and extrema:
adrcyc = AdrCyc(ap); % row-by-row cyclic permutations of ap
ecyc = zeros(1,per);
ecyc(1) = epr;
for k = 2:per
	ecyc(k) = FixedPt(AdrConv(adrcyc(k,:)), phi, p);
end

% Calculating all extrema:
eall = ecyc;
K = ceil(pi/abs(angle(phi(1))))+1; % needed iteration number for 180 deg
% Iterating to the left:
for k = 1:K
	eall = [ eall, p(1)+(phi(1)^k)*(ecyc-p(1)) ];
end
K = ceil(pi/abs(angle(phi(2))))+1;  % needed iteration number for 180 deg
for k = 1:K
	eall = [ eall, p(2)+(phi(2)^k)*(ecyc-p(2)) ];
end
% Keeping the relevant points:
chind = convhull(real(eall), imag(eall));
eall = eall(chind);
