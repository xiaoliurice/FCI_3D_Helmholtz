%% main script for the matrix-free preconditioner of the Helmholtz equation

% parametrization of the contour
np   = 8;    % num of poles
rad1 = 1;    % radius in x
rad2 = 5;    % radius in y
theta = pi/np.*(1:2:np*2-1);  % angles in radi

FCI = struct('np',np,'shf',[],'wts',[],'nim',8,'im',20,'tol',1e-4);
FCI.shf = complex( rad1*cos(theta), rad2*sin(theta) );             % points on the ellipse
FCI.wts = complex( rad2*cos(theta), rad1*sin(theta) )./FCI.shf/np; % weights
FCI.shf = FCI.shf - max(real(FCI.shf));                            % shift the center

MAT = struct('N',[],'kh0',[],'kh',[],'ab',[],'z',1); 

for rnd = 1
    N = ones(1,3)*rnd*40; % grid size
    MAT.N   = N;
    
    % absorbing layers
    nab = 8;    % number of points
    rab = 4e-4; % decay of amplitude
    cab = 1j*log(1/rab)^2;
    l1 = taper(N(1),nab,nab)*cab;
    l2 = taper(N(2),nab,nab)*cab;
    l3 = taper(N(3),nab,nab)*cab;
    MAT.ab  = zeros(N);
    MAT.ab(:) = max( kron(ones(N(2)*N(3),1), l1), kron(l3, ones(N(2)*N(1),1)) );
    MAT.ab(:) = max( MAT.ab(:), kron( ones(N(3),1), kron(l2, ones(N(1),1)) ) );
    
    % wavenumber
    MAT.kh0 = 2*pi/5.25;      % 2pi/number of points per wavelength
    MAT.kh  = ones(N)*MAT.kh0;
%    MAT.kh(ceil(N(1)/3):ceil(N(1)*2/3),ceil(N(2)/3):ceil(N(2)*2/3),ceil(N(3)/3):ceil(N(3)*2/3)) = MAT.kh0/2;
    
    
    % right hand side
    f = zeros(N);
    f(ceil(N(1)/2),ceil(N(2)/2),ceil(N(3)/2)) = 1;
    f = reshape(f,[prod(N),1]);
    nrm0 = norm(f);
    tol = 1e-6;
    
    % iterative refinement on FCI solution
    u   = zeros(size(f));
    nmv = 0;
    for i = 1:60
        [v,nm] = fcisol(f,MAT,FCI);
        nmv = nmv + nm;
        u = u + v;
        f = f - helmop(v,MAT);
        
        rres = norm(f)/nrm0;
        fprintf('it: %d, nm: %d, rres: %.2e\n\n',i,nm,rres);
        if (rres < tol)
            break;
        end
    end
    u = reshape(u,N);
    %imagesc(real(u(:,:,ceil(N(3)/2))));
end