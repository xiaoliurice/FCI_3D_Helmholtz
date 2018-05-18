function [sol, nmv, relres] = ric_poly(rhs,z,MAT,nblock, tol,niter, d)
% richardson iteration
% sol = sol + d*res

n = length(rhs);

if(nblock == 1)
    sol = zeros(size(rhs));
else
    rhs = [zeros(size(rhs));rhs];
    sol = zeros(2*n, 1);
end

res  = rhs;
nrm0 = norm(rhs);
for it = 1:niter
    sol = sol + d*res;
    res = rhs - helmop(sol,MAT) + z*sol;
    relres = norm(res) / nrm0;
    if(relres < tol)
        break;
    end
end


if(nblock == 2)
    sol = sol(n+1:end);
end
nmv = it;
end