%% main script for the matrix-free preconditioner of the Helmholtz equation
visual = false;  % whether to visualize
sparse = false;  % sparse or dense system

% refine
refine = 1;
if(sparse)
    refine = 4; % refine by 4 times
end

% a reference set of size, imaginary shift, and sampling rate
sz0  = 40;
sep0 = 0.4;
ppw0 = 2.25;

for sz = [160,240,320]
%for sz = [20,40,60,80]
    %% setup matrix
    N   = ones(1,3)*sz*refine; % grid size
    ppw = ppw0*refine;         % sampling rate
    MAT = mat_setup(N,pi/ppw,pi*2/ppw,sparse) % form the matrix
    if(visual)
        visualize(1./sqrt(MAT.M), [ceil(N(1)*0.3), ceil(N(2)*0.7), ceil(N(3)*0.3)],1);
    end
    
    %% solver
    sep = sep0*(sz0/sz)*refine; % (IMPORTANT) distance from real line
    tol = 1e-4;                 % tolerance
    np = 6;                     % num of poles
    for nblock = 1:2
        if(nblock == 1)
            np = 6;
        else
            np = 1;
        end
        for method = 1:2
            
            FCI = fci_setup(np,sep,0.5,nblock,method,MAT)
            
            % FCI iterative solution
            % right hand side
            f = zeros([prod(N),1]);
            f( ceil(N(1)/2) + N(1)*( floor(N(2)/2) + N(2)*floor(N(3)/2) ) ) = 1;
            nrm0 = norm(f);
            
            % iterative refinement on FCI solution
            u   = zeros(size(f));
            nmv  = 0;
            otimer = tic;
            for i = 1:1000
                [v,nmv1] = fcisol(f,MAT,FCI);
                nmv = nmv + nmv1 + 1;
                u = u + v;
                f = f - helmop(v,MAT);
                rres = norm(f)/nrm0;
                fprintf('it: %d, nmatvec: %d, relres: %.2e\n\n',i,nmv,rres);
                if (rres < tol); break; end
            end
            t = toc(otimer); fprintf('Total time: %.2f\n\n\n',t);
        end
        %sep = sep/sqrt(MAT.rho(1));
    end
    %%  visualization
    u = reshape(u(end-prod(N)+1:end),N);
    if(visual)
        if(sparse)
            % visualize
            fig2 = visualize(real(u),[ceil(N(1)*0.3), ceil(N(2)*0.7), ceil(N(3)*0.3)], 1);
        else
            % band-limited interpolation
            ninter = 3;
            tmp = zeros(N*ninter);
            m1 = floor(N(1)/2);
            m2 = floor(N(2)/2);
            m3 = floor(N(3)/2);
            tmp([1:m1, end-(N(1)-m1)+1:end],[1:m2, end-(N(2)-m2)+1:end],[1:m3, end-(N(3)-m3)+1:end]) = fftn(u);
            tmp = ifftn(tmp);
            fig2 = visualize(real(tmp),[ceil(N(1)*ninter*0.3), ceil(N(2)*ninter*0.7), ceil(N(3)*ninter*0.3)], ninter);
        end
    end
    
    %% unpreconditioned GMRES
%     otimer = tic;
%     [u,~,relres,iter] = gmres(@(x)helmop(x,MAT),f,FCI.im,tol,ceil(nmv/FCI.im) );
%     t = toc(otimer);
%     fprintf('GMRES: nmatvec: %d, relres: %.2e\n',(iter(1)-1)*FCI.im + iter(2),relres);
%     fprintf('Total time: %.2f\n\n',t);
end
