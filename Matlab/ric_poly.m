function [sol, nmv, relres] = ric_poly(rhs,z,MAT, tol,niter,d)
% richardson iteration
% sol = sol + d*res

res  = rhs;
sol  = zeros(size(rhs));
nrm0 = norm(rhs);
for it = 1:niter
    sol = sol + d*res;
    res = rhs - helmop(sol,MAT) + z*sol;
    relres = norm(res) / nrm0;
    if(relres < tol)
        break;
    end
end

nmv = it;
end