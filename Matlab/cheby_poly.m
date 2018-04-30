function [u, it] = cheby_poly(f, MAT, z, tol, itmax)
% Chebyshev iteration for solving MAT*u = f
% based on Algorithm 3 of [Gutknecht, Rollin, Parallel Computing 2002]

% the spectrum of MAT is the line from [-z,rho-z]
% contained in an ellipse with foci [a-c,a+c]
a = MAT.rho/2 - z;
c = MAT.rho/2;

r  = f./MAT.M(:);
u  = zeros(size(f));
du = zeros(size(f));

nrm0 = norm(r); % initial residual norm
relres = 1;     % current relative residual

it = 0;
beta = 0;
while(relres >= tol && it < itmax)
    gamma = -(a+beta);
    du = (-r + beta*du) / gamma;
    u  = u + du;
    r  = (f-helmop(u,z,MAT))./MAT.M(:); % new residual
    relres = norm(r)/nrm0;
    it = it + 1;
    if(it < 2)
        beta = -c*c/(2*a);
    else
        beta = (c/2)^2/gamma;
    end
end

end