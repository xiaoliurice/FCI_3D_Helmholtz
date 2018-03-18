%% main script for the matrix-free preconditioner of the Helmholtz equation
sz = [40,80,160];
visual = true;


for run = 2
    z0 = 1;
    N = ones(1,3)*sz(run);
    % set up the matrix: DST*S*DST - (z0+1j*D)*M
    % S, D, M are diagonal
    MAT = struct('N',N,'S',[],'M',[],'D',[], 'atmp',[],'btmp',[]);
    
    % eigenvalues of Laplaician in DST-I basis
    
    MAT.S = kron(ones(1,N(2)*N(3)),(pi/(N(1)+1) * (1:N(1))).^2) ...
        +kron( (pi/(N(3)+1) * (1:N(3))).^2, ones(1,N(1)*N(2)) )...
        +kron(ones(1,N(3)), kron( (pi/(N(2)+1) * (1:N(2))).^2, ones(1,N(1))));
    
    % squared wavenumber
    khmax = 2*pi/2.22;
    khmin = 2*pi/2.6;
    MAT.M  = ones(N)*khmax^2;
    MAT.M(ceil(N(1)/2):ceil(N(1)*2/3),ceil(N(2)/2):ceil(N(2)*2/3),ceil(N(3)/2):ceil(N(3)*2/3)) = khmin^2;
    MAT.M = smooth3(MAT.M,'gaussian',9,1);
    if(visual)
        fig1 = visualize(MAT.M, [ceil(N(1)*0.5), ceil(N(2)*0.5), ceil(N(3)*0.5)],1);
    end
    
    %absorbing layers
    nab = 10;    % number of points
    rab = 1e-2;  % decay of amplitude
    cab = abs(log(rab))/khmin*2;
    l1 = taper(N(1),nab,nab)*cab;
    l2 = taper(N(2),nab,nab)*cab;
    l3 = taper(N(3),nab,nab)*cab;
    MAT.D  = zeros(N);
    MAT.D(:) = max( kron(ones(N(2)*N(3),1), l1), kron(l3, ones(N(2)*N(1),1)) );
    MAT.D(:) = max( MAT.D(:), kron( ones(N(3),1), kron(l2, ones(N(1),1)) ) );
    clear nab rab cab l1 l2 l3
    
    MAT
    dmax = max(MAT.D(:))
    
    %% set up the solver
    np  = 2;       % num of poles
    rad = 0.2/2^(run-1);     % radius for real part
    sep = 0.2/2^(run-1);     % distance for imaginary part
    
    FCI = struct('np',np,'shf',[],'wts',[],'nim',1,'im',20,'tol',[0.3,0.001], 'dt',[],'num',[]);
    if(np==1)
        % shifted Laplacian
        FCI.shf = sep*1j;
        FCI.wts = 1;
    else
        theta = pi/np.*(1:2:np*2-1);           % angles in radi
        rad2 = (sep + dmax) / sin(theta(1));   % the first and last pole: (0,sep), (0,-sep-dmax*2)
        rad = rad2;                            % make a circle
        FCI.shf = complex( rad*cos(theta), rad2*sin(theta) );                  % points on the ellips
        FCI.shf = FCI.shf - complex(rad*cos(theta(1)), dmax);                  % shift the center
        FCI.wts = complex( rad2*cos(theta), rad*sin(theta) )./FCI.shf/np;      % weights
        fprintf('radius: %.2e, %.2e\n',rad,rad2);
        clear theta rad rad2 sep
    end
    FCI.num = zeros(size(FCI.wts));
    FCI.dt  = zeros(size(FCI.wts));
    rates   = zeros(size(FCI.wts));
    for p = 1:np
        z = z0+FCI.shf(p);
        if( imag(z) < 0 )
            z = z + 1j*dmax;
        end
        [FCI.dt(p), rates(p)] = convrate_exp(MAT.S(end)/khmin^2, z, 4);
        FCI.num(p) = ceil(log(FCI.tol(1)*abs(FCI.wts(1)/FCI.wts(p)))/log(rates(p)));
    end
    FCI.num = max(FCI.num, 1);
    FCI.nim = ceil(sum(FCI.num)*4/FCI.im );
    
    FCI
    rates
    
    %% iterative solution
    %right hand side
    f = zeros(N);
    f(ceil(N(1)/2),ceil(N(2)/2),ceil(N(3)/2)) = 1;
    f = reshape(f,[prod(N),1]);
    nrm0 = norm(f);
    tol = 1e-6;
    
    % iterative refinement on FCI solution
    u   = zeros(size(f));
    nmv = 0;
    otimer = tic;
    for i = 1:60
        [v,nmv1] = fcisol(f,z0,MAT,FCI);
        nmv = nmv + nmv1 + 1;
        u = u + v;
        f = f - helmop(v,z0,MAT);
        
        rres = norm(f)/nrm0;
        fprintf('it: %d, nmatvec: %d, rres: %.2e\n\n',i,nmv,rres);
        if (rres < tol)
            break;
        end
    end
    t = toc(otimer);
    fprintf('time: %.2f, nmatvec %d\n\n\n', t, nmv);
    u = reshape(u,N);
    
    %% band-limited interpolation and visualization
    if(visual)
        m = 3;             % refinement ratio
        tmp = zeros(N*m);
        tmp(1:N(1),1:N(2),1:N(3)) = dst1fft(u);
        tmp = dst1fft(tmp);
        fig2 = visualize(real(tmp), [ceil(N(1)*m*0.5), ceil(N(2)*m*0.5), ceil(N(3)*m*0.5)], m);
    end
end