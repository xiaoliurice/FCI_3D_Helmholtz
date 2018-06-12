function [u,nmv] = fcisol(f,MAT,FCI)
% one step of contour integration solution
% because of the last step of GMRES, outside this function an iterative refinement suffices

nrm0 = norm(f);
n = length(f);
u = zeros(size(f));

if(FCI.nblock==2)
    f = [zeros(size(f));f];
end

ftimer = tic;
% outer problem
nmv = 0;
nmvp = 0;
for p = 1:FCI.np
    z = FCI.shf(p);
    tol = FCI.tol(1) * sqrt(abs(FCI.wts(1)/FCI.wts(p)));
    
    if(FCI.method==1)
        [v,nmv1, relres] =  exp_poly(f,z,MAT, tol, FCI.num(p), complex(FCI.bet(1),imag(z)),FCI.d(p),FCI.q(p));
        fprintf('|%.1e,%d',relres*abs(FCI.wts(p))/abs(FCI.wts(1)),nmv1);
    else
        [v,nmv1, relres] = cheby_poly(f,z,MAT, tol, FCI.num(p)*FCI.q(p));
        fprintf('|%.1e,%d',relres*abs(FCI.wts(p))/abs(FCI.wts(1)),nmv1);       
        w = f - (helmop(v,MAT) - z*v);
        nrm = norm(w);
        if(nrm > nrm0*tol)
            [w,~,relres,iter] = gmres(@(x)helmop(x,MAT)-z*x, w, FCI.im, nrm0*tol/nrm, FCI.num(p)*FCI.q(p));
            v = v + w;
            nmv1 = nmv1 + (iter(1)-1)*FCI.im + iter(2);
            fprintf(';%.1e,%d', nrm*relres/nrm0 * abs(FCI.wts(p))/abs(FCI.wts(1)), nmv1);
        end
    end
    nmvp = nmvp + nmv1;
    
    % project into the original space
    if(FCI.nblock==2)
        v = v(n+1:n*2);
    end
    u = u + FCI.wts(p)*v;
end
nmv = nmv + nmvp;

v = helmop(u,MAT);
c = (v'*f(end-n+1:end))/(v'*v); % step size
u = c*u;
v = f(end-n+1:end)-c*v; % new residual
fprintf('|outer step (%.2f,%.2f)',real(c),imag(c));

% inner problem
% GMRES on current residual
[v,~,relres,iter] = gmres(@(x)helmop(x,MAT),v,FCI.im,FCI.tol(2),ceil(nmvp*2/FCI.im));
u = u + v;
nmv1 = (iter(1)-1)*FCI.im + iter(2) + 1;
nmv = nmv + nmv1;
fprintf('|inner %.2e, %d',relres, nmv1);

t=toc(ftimer); fprintf('|time %.2fs\n',t);
end
