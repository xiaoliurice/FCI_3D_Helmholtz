%% main script for the matrix-free preconditioner of the Helmholtz equation
visual = true; % whether to visualize
sparse = false; % sparse or dense system

% a reference set of size, imaginary shift, and sampling rate
sz0  = 40;
sep0 = 0.2;
ppw0 = 2.22;

for sz = sz0*[1,2]
    %% setup the matrix
    N   = ones(1,3)*sz; % grid size
    ppw = ppw0;         % sampling rate
    if(sparse)
        % refine by 4 times for 7-pt stencil
        N   = N*4;
        ppw = ppw*4;
    end
    MAT = mat_setup(N,pi/ppw,pi*2/ppw,sparse) % form the matrix

    if(visual)
        visualize(MAT.M, [ceil(N(1)*0.3), ceil(N(2)*0.7), ceil(N(3)*0.3)],1);
    end
    
    % set up the solver
    np  = 1;                        % num of poles
    sep = sep0*(sz0/sz); % (IMPORTANT) distance from real line
    FCI = fci_setup(np,sep,0.2,MAT)
    
    %% iterative solution
    %right hand side
    f = zeros([prod(N),1]);
    f( ceil(N(1)/2) + N(1)*( floor(N(2)/2) + N(2)*floor(N(3)/2) ) ) = 1;
    nrm0 = norm(f);
    
    % iterative refinement on FCI solution
    u   = zeros(size(f));
    nmv = 0;
    otimer = tic;
    for i = 1:60
        [v,nmv1] = fcisol(f,MAT,FCI);
        nmv = nmv + nmv1 + 1;
        u = u + v;
        f = f - helmop(v,1,MAT);
        rres = norm(f)/nrm0;
        fprintf('it: %d, nmatvec: %d, relres: %.2e\n\n',i,nmv,rres);
        if (rres < 1e-6); break; end
    end
    t = toc(otimer); fprintf('Total time: %.2f\n\n\n',t);

    
    %%  visualization
    u = reshape(u,N);
    if(visual)
        if(sparse)
            % visualize
            fig2 = visualize(real(u),[ceil(N(1)*0.3), ceil(N(2)*0.7), ceil(N(3)*0.3)], m);
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
end