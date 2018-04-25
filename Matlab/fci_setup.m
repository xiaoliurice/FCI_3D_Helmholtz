function [FCI] = fci_setup(np,sep,asp,MAT)
% set up the solver
% np   = num of poles
% sep  = min separation from the real axis
% asp  = aspact ratio of the contour


% 'shf' shifts, 'wts' weights
% 'nim' [num of restarts, Krylov subspace dimension]
% 'tol' [tol of outer problems, tol of inner problem]
% 'dt' step size, 'num' num of fixed pnt iters, 'q' order
FCI = struct('np',np,'shf',[],'wts',[],'nim',[1,20],'tol',[2e-1,2e-1],'dt',[],'num',[],'q',4);

dmax = max(MAT.D(:))

% quadrature points and weights
if(np==1)
    % single imaginary shift
    FCI.shf = sep*1j;
    FCI.wts = 1;
else
    phi = pi/np.*(1:2:np*2-1);   % angles in radi
    dy  = dmax*0.7;              % shift  in y
    ry  = (sep+dy)/sin(phi(1));  % radius in y
    rx  = ry*asp;                % radius in x
    dx  = rx*cos(phi(1));        % shift  in x
    
    FCI.shf = complex(rx*cos(phi),ry*sin(phi)) - complex(dx,dy);  % points on the ellipse
    FCI.wts = complex(ry*cos(phi),rx*sin(phi))./FCI.shf/np;       % weights
    fprintf('radius: %f, %f, center: %f,%f\n',rx,ry,-dx,-dy);
end

% set up fixed point iterations
FCI.num = zeros(size(FCI.wts));
FCI.dt  = zeros(size(FCI.wts));
rates   = zeros(size(FCI.wts));
for p = 1:np
    z = 1+FCI.shf(p);
    if( imag(z) < 0 )
        % distance in y from the nearest eigenvalue
        z = z + 1j*dmax;
    end
    [FCI.dt(p), rates(p)] = exp_rate(MAT.rho, z,FCI.q);
    FCI.num(p) = ceil( log(FCI.tol(1) * abs(FCI.wts(1)/FCI.wts(p)).^0.5 )/log(rates(p)));
end
FCI.num = max(FCI.num, 1);
FCI.nim(1) = ceil(sum(FCI.num)*FCI.q/FCI.nim(2) );

rates
end