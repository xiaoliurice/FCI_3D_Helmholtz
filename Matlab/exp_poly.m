function [sol, nmv, relres] = exp_poly(rhs,z,MAT, tol,niter, z0,d,q)
% solution of shifted Helmholtz via p-th order taylor expansion of exp func

wd  = -1j*d;
wdz = wd*(z-z0);
wv  = ones(1,q+1);
for j=1:q
    wv(j+1) = wv(j) * (wdz/j);
end
wsum = sum(wv);

sol = zeros(size(rhs));
nrm0 = norm(rhs);
nmv = 0;
for i = 1:niter
    k = sol;
    for j=1:q
        k = (wd/j) * (helmop(k,MAT) - z0*k - wv(j)*rhs);
        nmv = nmv + 1;
        if(j==1)
            relres = norm(k-wdz*sol)/nrm0/abs(wd);
            if(relres<tol)
                break;
            end
        end
        sol  = sol + k;
    end
    if(relres<tol)
        break;
    end
    sol  =  sol/wsum;
end
end