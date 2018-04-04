function [MAT] = mat_setup(N,khmin,khmax)
    % set up the matrix: DST*S*DST - (z0+1j*D)*M
    % S, D, M are diagonal
    MAT = struct('N',N,'S',[],'M',[],'D',[], 'atmp',[],'btmp',[]);
    
    % eigenvalues of Laplaician in DST-I basis
    MAT.S = kron(ones(1,N(2)*N(3)),(pi/(N(1)+1) * (1:N(1))).^2) ...
        +kron( (pi/(N(3)+1) * (1:N(3))).^2, ones(1,N(1)*N(2)) )...
        +kron(ones(1,N(3)), kron( (pi/(N(2)+1) * (1:N(2))).^2, ones(1,N(1))));
    
    % squared wavenumber
    MAT.M  = ones(N)*khmax^2;
    MAT.M(ceil(N(1)/2):ceil(N(1)*2/3),ceil(N(2)/2):ceil(N(2)*2/3),ceil(N(3)/2):ceil(N(3)*2/3)) = khmin^2;
    MAT.M = smooth3(MAT.M,'gaussian',9,1);
    
    %absorbing layers
    nab = 10;    % number of points
    rab = 1e-2;  % decay of amplitude
    cab = abs(log(rab))/khmin*2;
    l1 = taper(N(1),nab,nab)*cab;
    l2 = taper(N(2),nab,nab)*cab;
    l3 = taper(N(3),nab,nab)*cab;
    MAT.D  = zeros(N);
%     MAT.D(:) = max( kron(ones(N(2)*N(3),1), l1), kron(l3, ones(N(2)*N(1),1)) );
%     MAT.D(:) = max( MAT.D(:), kron( ones(N(3),1), kron(l2, ones(N(1),1)) ) );
end