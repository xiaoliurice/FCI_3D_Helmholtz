%% script for checking the polynomial solution of shifted problems
visual = true;
sparse = false;

sz = 40;

N = ones(1,3)*sz; % grid size
ppw = 2.22;       % sampling rate
if(sparse)
    % refine by 4 times for 7-pt stencil
    N   = N*4;
    ppw = ppw*4;
end
MAT = mat_setup(N,pi/ppw,pi*2/ppw,sparse) % form the matrix
    
% right hand side
f = rand([prod(N),1]);
nrm0 = norm(f);
tol = 1e-3;

%% select step size and perform a trial run

% compute the square root for splitting
for z = [1+1j] %, 1+0.5j, 1+0.25j, 1+ 0.125j, 1+0.0625j]
    ss = sign(imag(z));
    tmp = sqrt(z + 1j*MAT.D(:));
    MAT.atmp = real(tmp).^2;
    MAT.btmp = (ss*imag(tmp)) ./ real(tmp);
    
    for q = 4
        [dt, rate] = exp_rate(MAT.rho,z,q);
        num = ceil( log(tol)/log(rate) );
        disp('--order, step size, rate per matvec, rate, nmatvec')
        fprintf('%d, %.3f, %.3f, %.2e, %d\n',...
                q,dt, rate.^(1/q), rate, num*q);
        
        % true run
        [u,nmv] = exp_poly(f,MAT,ss,num,dt,q);
        fprintf('relres %.2e\n',norm(f-helmop(u,z,MAT)) / nrm0);
    end
end