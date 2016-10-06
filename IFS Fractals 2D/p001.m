% Tests if min_c max_k |p_k - c| = max_k |p_k - ave(p)|.

clear all;

n = 7;

% Gen main pts
p = 100*(rand(n,1) + rand(n,1)*i);

% Gen grid
grd = zeros(10000,2);
for j = 1:100
	for k = 1:100
		grd(100*(j-1)+k,1) = j+k*i;
	end
end

% Compute values
for j = 1:10000
	c = grd(j,1);
	grd(j,2) = p001f(p,c);
end

% Find minimum
minf = min(grd(:,2));
disp(grd(find(grd(:,2) == minf)',1));
disp(['f(Min) = ', num2str(minf)])

% Value at the center of mass
disp(['Com = ', num2str(mean(p))])
disp(['f(CoM) = ', num2str(p001f(p,mean(p)))])

% Value at harmonic center
fpj = p001f(p,p(1));
for j = 2:n
	fpj = [ fpj, p001f(p,p(j))];
end
ifpj = 1./fpj;
hc = sum(ifpj*p)/sum(ifpj);
disp(['HC = ', num2str(hc)]);
disp(['f(HC) = ', num2str(p001f(p,hc))])

% fminsearch
[x,fval] = fminsearch(@(c)max(abs(p-c)),5+5*i);
disp(['c = ', num2str(x)])
disp(['f(c) = ', num2str(fval)])

