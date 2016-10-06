% Principal extremal point and address for C-type bifractals (theoretically).
% Need: phi(1)<0, phi(2)>0
% epr =		principal extremal point
% adrpr =	address of the pr ext pt
% ecyc =	cyclical extremal pts (incl epr)
% adrcyc =	addresses of oth ext pts
% eoth =	all ext pts (sorted angularly)
% adroth =	addresses of all ext pts

function [ epr, ecyc, eall, adrpr, adrcyc, adrall ] = ExtPts2(phi, p, L)

adr = zeros(2,L);

% Calculating epr and adr(epr)
r = -angle(phi(2))/angle(phi(1));
if (r > 1)
	mp = [ 1 2 ];
else
	r=1/r;
	mp = [ 2 1 ];	
end
n = 0;
s = 0;
for l = 1:L
	%n = max(0, round(r*l-s));
	n = round(r*l-s);
	s = s+n;
	adr(1,2*l-1) = mp(2);
	adr(2,2*l-1) = 1;
	adr(1,2*l) = mp(1);
	adr(2,2*l) = n;
	
end
epr = AdrB(p(1), adr, phi, p);
adrpr = AdrConv(adr);

% Calculating the cyclical extrema and addresses
Lv = length(adrpr);
[ per, ap ] = SeqPeriod(adrpr);
adrcyc0 = AdrCyc(ap);
ecyc = epr;
adrcyc = adrpr;
for k = 2:per
	adrcyc = [ adrcyc; zeros(1,Lv) ];
	for j = 1:Lv
		col = mod(j,per);
		if (col == 0) col=per; end
		adrcyc(k,j) = adrcyc0(k, col);
	end
	ecyc = [ ecyc, AdrB(p(1), [adrcyc(k,:); ones(1,Lv)], phi, p) ];
end

% Calculating all extrema

%arg1 = max(angle(ecyc-p(1)))-abs(angle(phi(1)));
%relind = find(angle(ecyc-p(1)) > arg1);
%erel = ecyc(relind);
erel = ecyc;
eall = erel;
nrel = length(erel);
% adrrel = adrcyc(relind,:);
adrrel = adrcyc;
adrrel0 = adrrel;
adrall = adrrel;
K = ceil(pi/abs(angle(phi(1))))+1;
for k = 1:K
	eall = [ eall, p(1)+(phi(1)^k)*(erel-p(1)) ];
	adrrel = [ 1*ones(nrel,k), adrrel0(:,1:(Lv-k)) ];
	adrall = [ adrall; adrrel ];
end

%arg2 = min(angle(ecyc-p(2)))+abs(angle(phi(2)));
%relind = find(angle(ecyc-p(2)) < arg2);
%erel = ecyc(relind);
erel = ecyc;
nrel = length(erel);
%adrrel = adrcyc(relind,:);
adrrel = adrcyc;
adrrel0 = adrrel;
K = ceil(pi/abs(angle(phi(2))))+1;
for k = 1:K
	eall = [ eall, p(2)+(phi(2)^k)*(erel-p(2)) ];
	adrrel = [ 2*ones(nrel,k), adrrel0(:,1:(Lv-k)) ];
	adrall = [ adrall; adrrel ];
end

% Keeping the relevant points
[ eall, ind ] = ConvHullTol(eall, 0.001);
adrall = adrall(ind,:);
