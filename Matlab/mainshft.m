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
% the following line removes absorbing boundaries
%MAT.D = MAT.D*0;

% right hand side
f = rand([prod(N),1]);
%f = zeros([prod(N),1]);
%f( ceil(N(1)/2) + N(1)*( floor(N(2)/2) + N(2)*floor(N(3)/2) ) ) = 1;
nrm0 = norm(f);
tol = 1e-5;

%% select step size and perform a trial run
% compute the square root for splitting
for z = [1+1j, 1+0.5j, 1+0.25j, 1+ 0.125j, 1+0.0625j]
    for q = 5
        [dt, rate] = exp_rate(MAT.rho,z,q);
        num = ceil( log(tol)/log(rate) );
        disp('--order, step size, rate per matvec, rate, nmatvec')
        fprintf('%d, %.3f, %.3f, %.2e, %d\n',...
                q,dt, rate.^(1/q), rate, num*q);
        
        % true runs
        [u,nmv] = exp_poly(f,MAT,z,num,dt,q);
        fprintf('Taylor of exp: relres %.2e, matvecs %d\n',norm(f-helmop(u,z,MAT)) / nrm0, nmv);
        
        [u,nmv] = cheby_poly(f,MAT,z,tol,num*q);
        fprintf('Chebyshev: relres %.2e, matvecs %d\n',norm(f-helmop(u,z,MAT)) / nrm0, nmv);
        
        [u,nmv] = cheby_poly2(f,MAT,z,tol,num*q);
        fprintf('Chebyshev on 2x2 system: relres %.2e, matvecs %d\n',norm(f-helmop(u,z,MAT)) / nrm0, nmv);
    end
end

%%
u = reshape(u,N);
if(visual)
    if(sparse)
        % visualize
        fig2 = visualize(real(u),[ceil(N(1)*0.5), ceil(N(2)*0.5), ceil(N(3)*0.5)], m);
    else
        % band-limited interpolation
        ninter = 3;
        tmp = zeros(N*ninter);
        m1 = floor(N(1)/2);
        m2 = floor(N(2)/2);
        m3 = floor(N(3)/2);
        tmp([1:m1, end-(N(1)-m1)+1:end],[1:m2, end-(N(2)-m2)+1:end],[1:m3, end-(N(3)-m3)+1:end]) = fftn(u);
        tmp = ifftn(tmp);
        fig2 = visualize(real(tmp),[ceil(N(1)*ninter*0.5), ceil(N(2)*ninter*0.5), ceil(N(3)*ninter*0.5)], ninter);
    end
end