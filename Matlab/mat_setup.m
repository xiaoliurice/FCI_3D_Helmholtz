function [MAT] = mat_setup(N,khmin,khmax,sparse)
    % set up the matrix: stiffness - (z0+1j*D)*M
    % N = grid size
    % khmin = 2pi / maximum sampling rate
    % khmax = 2pi / minimum sampling rate
    % sparse = 1 for 7-point stencil  (sampling rate > 8);
    %        = 0 for Fourier spectral (sampling rate > 2).
    
    % 'S' stiffness, 'M' mass, 'D' boundary, 'rho' rho(inv(M)S)
    MAT = struct('N',N, 'S',[], 'M',[], 'D',[], 'sparse',sparse, 'rho', 0);
    
    % mass matrix: squared wavenumber
    MAT.M  = ones(N)*khmax^2;
    MAT.M(ceil(N(1)*0.2):ceil(N(1)*0.4),ceil(N(2)*0.2):ceil(N(2)*0.4),ceil(N(3)*0.2):ceil(N(3)*0.4)) = khmin^2;
    MAT.M(ceil(N(1)*0.2):ceil(N(1)*0.4),ceil(N(2)*0.2):ceil(N(2)*0.4),ceil(N(3)*0.6):ceil(N(3)*0.8)) = khmin^2;
    MAT.M(ceil(N(1)*0.2):ceil(N(1)*0.4),ceil(N(2)*0.6):ceil(N(2)*0.8),ceil(N(3)*0.2):ceil(N(3)*0.4)) = khmin^2;
    MAT.M(ceil(N(1)*0.2):ceil(N(1)*0.4),ceil(N(2)*0.6):ceil(N(2)*0.8),ceil(N(3)*0.6):ceil(N(3)*0.8)) = khmin^2;
    MAT.M(ceil(N(1)*0.6):ceil(N(1)*0.8),ceil(N(2)*0.2):ceil(N(2)*0.4),ceil(N(3)*0.2):ceil(N(3)*0.4)) = khmin^2;
    MAT.M(ceil(N(1)*0.6):ceil(N(1)*0.8),ceil(N(2)*0.2):ceil(N(2)*0.4),ceil(N(3)*0.6):ceil(N(3)*0.8)) = khmin^2;
    MAT.M(ceil(N(1)*0.6):ceil(N(1)*0.8),ceil(N(2)*0.6):ceil(N(2)*0.8),ceil(N(3)*0.2):ceil(N(3)*0.4)) = khmin^2;
    MAT.M(ceil(N(1)*0.6):ceil(N(1)*0.8),ceil(N(2)*0.6):ceil(N(2)*0.8),ceil(N(3)*0.6):ceil(N(3)*0.8)) = khmin^2;
    MAT.M = smooth3(MAT.M,'gaussian',9,1);
    
    % non-hermitian part from absorbing layers, set cab=0 to remove it
    nab = 8;     % number of points
    if(sparse)
        nab = 20;
    end
    rab = 1e-2;  % decay of amplitude
    cab = abs(log(rab))/khmax*2;
    l1 = taper(N(1),nab,nab)*cab;
    l2 = taper(N(2),nab,nab)*cab;
    l3 = taper(N(3),nab,nab)*cab;
    MAT.D  = zeros(N);
    MAT.D(:) = max( kron(ones(N(2)*N(3),1), l1), kron(l3, ones(N(2)*N(1),1)) );
    MAT.D(:) = max( MAT.D(:), kron( ones(N(3),1), kron(l2, ones(N(1),1)) ) );
    
    % stiffness matrix
    if(sparse)
        e1 = ones(N(1),1);
        e2 = ones(N(2),1);
        e3 = ones(N(3),1);
        T1 = spdiags([-e1 2*e1 -e1], -1:1, N(1), N(1));
        T2 = spdiags([-e2 2*e2 -e2], -1:1, N(2), N(2));
        T3 = spdiags([-e3 2*e3 -e3], -1:1, N(3), N(3));
        Id1 = spdiags(e1, 0, N(1), N(1));
        Id2 = spdiags(e2, 0, N(2), N(2));
        Id3 = spdiags(e3, 0, N(3), N(3));
        MAT.S  = kron(kron(T3,Id2),Id1) + kron(kron(Id3,T2),Id1) + kron(kron(Id3,Id2),T1);
        MAT.rho =  12 / khmin^2;
    else
        % eigenvalues of Laplaician in DFT basis
        l1 = min(0:N(1)-1,N(1):-1:1) * (2*pi/N(1));
        l2 = min(0:N(2)-1,N(2):-1:1) * (2*pi/N(2));
        l3 = min(0:N(3)-1,N(3):-1:1) * (2*pi/N(3));
        MAT.S = kron(ones(1,N(2)*N(3)),l1.^2) ...
            +kron( l3.^2, ones(1,N(1)*N(2)) )...
            +kron(ones(1,N(3)), kron( l2.^2, ones(1,N(1))));
        MAT.rho = max(MAT.S(:)) / khmin^2;
    end
end

function v = taper(N, l1, l2)
% cosine taper at both end
v = zeros(N,1);
for i = 1:l1
    v(i) = 1+cos( i/(l1+1)*pi );
end
v(1:l1) = v(1:l1)/sum(v(1:l1));
for i = 1:l2
    v(N+1-i) = 1+cos( i/(l2+1)*pi );
end
v(N+1-l2:N) = v(N+1-l2:N)/sum(v(N+1-l2:N));
end