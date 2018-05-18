%% script for checking the polynomial solution of shifted problems
visual = true;
sparse = false;

refine = 1;
if(sparse)
    refine = 4; % refine by 4 times
end

sz = 40;
ppw = 2.25*refine;         %sampling rate
N   = ones(1,3)*sz*refine; %grid size

MAT = mat_setup(N,pi/ppw,pi*2/ppw,sparse) % form the matrix
MAT.rho = MAT.rho*1.2;
% the following line removes absorbing boundaries
%MAT.D = MAT.D*0; MAT.rho(2) = 0;

% right hand side
f = rand([prod(N),1]);
%f = zeros([prod(N),1]);
%f( ceil(N(1)/2) + N(1)*( floor(N(2)/2) + N(2)*floor(N(3)/2) ) ) = 1;
nrm0 = norm(f);
tol = 1e-2;


%% select step size and perform a trial run
% compute the square root for splitting

% upper part
vz = 1j*[1,0.5,0.25,0.125]*refine;
% lower part
%vz = - 1j*([1,0.5,0.25,0.125]*refine + MAT.rho(2));

for nblock = 1:2
    fprintf('\nnblock=%d\n',nblock);
    for z = vz
        % richardson
        [d, rate] = ric_rate(z,MAT.rho,nblock);
        num = ceil( log(tol)/log(rate) );
        disp('--richardson: rate, nmatvec, relres, nmatvec')
        fprintf('%.4f, %d\n',rate,num);
        [u, nmv, relres] = ric_poly(f,z,MAT,nblock, tol,num*2, d);
        fprintf('  relres %.2e, matvecs %d\n', relres, nmv);
        
        % Taylor expansion of matrix exponential
        disp('--order, step size, rate per matvec, rate, nmatvec')
        for q = 1:6
            [z0, d, rate] = exp_rate(z,MAT.rho,q,nblock);
            num = ceil( log(tol)/log(rate) );    
            fprintf('%d, %.3f, %.4f, %.4f, %d\n',...
                    q, d, rate.^(1/q), rate, num*q);
            % true runs
            [u, nmv, relres] = exp_poly(f,z,MAT,nblock, tol,num*2, z0,d,q);
            fprintf('  relres %.2e, matvecs %d\n', relres, nmv);
        end
        
    end
end

%%
% u = reshape(u,N);
% if(visual)
%     if(sparse)
%         % visualize
%         fig2 = visualize(real(u),[ceil(N(1)*0.5), ceil(N(2)*0.5), ceil(N(3)*0.5)], m);
%     else
%         % band-limited interpolation
%         ninter = 3;
%         tmp = zeros(N*ninter);
%         m1 = floor(N(1)/2);
%         m2 = floor(N(2)/2);
%         m3 = floor(N(3)/2);
%         tmp([1:m1, end-(N(1)-m1)+1:end],[1:m2, end-(N(2)-m2)+1:end],[1:m3, end-(N(3)-m3)+1:end]) = fftn(u);
%         tmp = ifftn(tmp);
%         fig2 = visualize(real(tmp),[ceil(N(1)*ninter*0.5), ceil(N(2)*ninter*0.5), ceil(N(3)*ninter*0.5)], ninter);
%     end
% end