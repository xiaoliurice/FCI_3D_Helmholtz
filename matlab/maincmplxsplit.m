%%
% the matrix is DST*DD*DST - 1j*BD - z*M (or M0)
% DD, BD, M are diagonal
MAT = struct('N',[],'DD',[],'M0',[],'M',[],'BD',[],'zf',1,'zs',1);
N = ones(1,3)*100; % grid size
MAT.N   = N;

% eigenvalues of Laplaician in DST-I basis
MAT.DD = kron(ones(1,N(2)*N(3)),(pi/(N(1)+1) * (1:N(1))).^2) ...
    +kron( (pi/(N(3)+1) * (1:N(3))).^2, ones(1,N(1)*N(2)) )...
    +kron(ones(1,N(3)), kron( (pi/(N(2)+1) * (1:N(2))).^2, ones(1,N(1))));

% squared wavenumber
khmax = 2*pi/2.22;
khmin = 2*pi/2.22;
MAT.M0 = khmax^2;      % square of 2pi/number of points per wavelength
MAT.M  = ones(N)*MAT.M0;
MAT.M(ceil(N(1)/2):ceil(N(1)*2/3),ceil(N(2)/2):ceil(N(2)*2/3),ceil(N(3)/2):ceil(N(3)*2/3)) = khmin^2;

% absorbing layers
nab = 8;    % number of points
rab = 1e-4; % decay of amplitude
cab = abs(log(rab))/khmin;
l1 = taper(N(1),nab,nab)*cab;
l2 = taper(N(2),nab,nab)*cab;
l3 = taper(N(3),nab,nab)*cab;
MAT.BD  = zeros(N);
MAT.BD(:) = max( kron(ones(N(2)*N(3),1), l1), kron(l3, ones(N(2)*N(1),1)) );
MAT.BD(:) = max( MAT.BD(:), kron( ones(N(3),1), kron(l2, ones(N(1),1)) ) );
clear nab rab cab l1 l2 l3

% right hand side
f = zeros(N);
f(ceil(N(1)/2)-8,ceil(N(2)/2)-8,ceil(N(3)/2)-8) = 1;
f = reshape(f,[prod(N),1]);
nrm0 = norm(f);
tol = 1e-6;

MAT.zs = 1 - 1j;
MAT.zf = MAT.zs;
MAT

%% select step size and perform a trial run
for p = 3
    [dt, rate] = convrate_exp(MAT.DD(end)/khmin^2,MAT.zs,p);
    disp('--order, step size, rate per matvec, rate, rate^10')
    [p,dt, rate.^(1/p),rate,rate^10]
    
    % true run
    num = 10;
    [u,nmv] = cmplxsplit_exp(f,MAT,num,dt,p);
    fprintf('relres %f\n',norm(f-helmop(u,MAT)) / nrm0);
end


