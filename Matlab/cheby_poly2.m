function [u, it] = cheby_poly2(f, MAT, z, tol, itmax)
% another version using a block system

ss   = sign(imag(z));
w    = -1j*ss;
tmp  = sqrt(z + 1j*MAT.D(:));
atmp = real(tmp).^2;
btmp = ss*imag(tmp) ./ real(tmp);

% foci [a-c,a+c]
tmp = sqrt(z);
a = w + ss*imag(tmp)/real(tmp);
c = 1j*sqrt(MAT.rho/real(tmp)^2);

r2 = f./MAT.M(:)./atmp;
r  = zeros(size(f));

u2 = zeros(size(f));
u  = zeros(size(f));

du2= zeros(size(f));
du = zeros(size(f));

nrm0 = norm(r2); % initial residual norm
relres = 1;      % current relative residual

it = 0;
beta = 0;
while(relres >= tol && it < itmax)
    gamma = -(a+beta);
    
    du  = (-r  + beta*du ) / gamma;
    du2 = (-r2 + beta*du2) / gamma;
    
    u  = u  + du;
    u2 = u2 + du2;
    
    % residual (r,r2)' = (0,(atmp*M)\f) - sysop(u,u2) - w(u,u2)'
    
    [v,v2] = sysop(u,u2,f,atmp,btmp, MAT);    
    r  = v  - w*u;
    r2 = v2 - w*u2;
    
    relres = sqrt(norm(r)^2+norm(r2)^2) / nrm0;
    
    it = it + 1;
    if(it < 2)
        beta = -c*c/(2*a);
    else
        beta = c*c/(4*gamma);
    end
end

end

function [y1,y2] = sysop(x1,x2,f,atmp,btmp, MAT)
% 2x2 system operator
%  / -b     1 \ /x1\   / 0   \
%  \-M\S/a  -b/ \x2/ + \M\f/a/

y1 = x2 - btmp.*x1;
y2 = (f - stiffop(x1,MAT))./MAT.M(:)./atmp - btmp.*x2;
end