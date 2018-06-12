function [sol, nmv, relres] = cheby_poly(rhs,z,MAT, tol,niter)
% Chebyshev iteration for solving MAT*sol = rhs
% based on Algorithm 3 of [Gutknecht, Rollin, Parallel Computing 2002]

% the spectrum of MAT is the line from [-1-z,rho-1-z]
% contained in an ellipse with foci [a-c,a+c]

if(prod(MAT.N) == length(rhs))
    c = MAT.rho(1)/2;
    a = c - (1+z);
else
    c = sqrt(MAT.rho(1));
    a = -(1+z);
end

r = rhs;
sol  = zeros(size(r));
dsol = zeros(size(r));

nrm0 = norm(r); % initial residual norm
relres = 1;     % current relative residual

nmv = 0;
beta = 0;
while(relres >= tol && nmv < niter)
    gamma = -(a+beta);
    dsol = (-r + beta*dsol) / gamma;
    sol  = sol + dsol;
    r  = rhs - (helmsym(sol,MAT) - z*sol); % new residual
    relres = norm(r)/nrm0;
    nmv = nmv + 1;
    if(nmv < 2)
        beta = -c*c/(2*a);
    else
        beta = (c/2)^2/gamma;
    end
end
end

function [y] = helmsym(x,MAT)
% Helmholtz operator
    n = prod(MAT.N);
    if(n==length(x))
        y = stiffop(x,MAT) - MAT.M(:).*x;
    else
        % double sized one
        y = 1j*[x(n+1:end);-stiffop(x(1:n),MAT)-(1-MAT.M(:)).*x(1:n)] - x;
    end
end