function iii = ints2(gaps, pots, temp, scat)

%  function for evaluating the six gap equations
%  and four number equations
%  usage: result = ints(gaps, pots, temp, scat)
%  returns result which contains the integrals,
%  result(1)-result(6) contain gap integral equations
%  result(7)-result(10) contain number equation ones
%  gaps should contain the given gaps (6)
%  pots for chemical potentials (4)
%  temp for 1/kT
%  scat for scattering lengths times fermi wavenumbers


gaps = reshape(gaps, [length(gaps) 1]);
pots = reshape(pots, [length(pots) 1]);
scat = reshape(scat, [length(scat) 1]);

potentialmatrix = diag([pots(1) pots(2) 0 pots(3) 0 pots(4) 0]);

gapmatrix = zeros(7,7);
gapmatrix(1,3) = gaps(1);
gapmatrix(1,5) = gaps(2);
gapmatrix(1,7) = gaps(3);
gapmatrix(2,5) = gaps(4);
gapmatrix(2,7) = gaps(5);
gapmatrix(5,7) = gaps(6);
gapmatrix = gapmatrix+gapmatrix';

discn = 5e4; cutoff = 30;
epoints = linspace(0,cutoff,discn);
sqepoints = sqrt(epoints);
renorm = 1./(2*sqepoints);

xisbase = [1 1 0 1 0 1 0];
xis = diag(xisbase);

Ks = zeros(7,7,discn);
Es = zeros(7,discn);
tempE = zeros(7,7);
Us = zeros(7,7,discn);

for iter = [1:discn]
	Ks(:,:,iter) = epoints(iter)*xis+gapmatrix-potentialmatrix;
	[Us(:,:,iter) tempE] = eig(Ks(:,:,iter));
	Es(:,iter) = diag(tempE);
end

fermis = (exp(temp*Es)+1).^(-1);
Gs = kron(ones(1,discn),xisbase')+((-2*xis+eye(7))*fermis);

%  plot(epoints,Gs)

iii = zeros(10,1);



iii(7)=cutoff*sum(sqepoints.*sum(abs(reshape(Us(:,1,:),[7 discn])).*(1-Gs)))/discn;
%  iii(8)=cutoff*sum(sqepoints.*sum(abs(reshape(Us(:,2,:),[7 discn])).*(1-Gs)))/discn;
%  iii(9)=cutoff*sum(sqepoints.*sum(abs(reshape(Us(:,4,:),[7 discn])).*(1-Gs)))/discn;
%  iii(10)=cutoff*sum(sqepoints.*sum(abs(reshape(Us(:,6,:),[7 discn])).*(1-Gs)))/discn;



Tcoef = temp

gridnum = discn
mu = pots(1)

points = linspace(0,cutoff,gridnum);
sqpoints = sqrt(points);

frm=(exp(Tcoef*(points-mu))+1).^(-1);

iii(9) = cutoff*sum(sqepoints.*frm)/gridnum;

