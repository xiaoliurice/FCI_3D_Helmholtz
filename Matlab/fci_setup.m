function [FCI] = fci_setup(np,sep,asp,nblock,method,MAT)
% set up the solver
% np   = num of poles
% sep  = min separation from the real axis
% asp  = aspact ratio of the contour
% nblock = 1 original; 2 extended
% method = 1 exp(); 2 Chebyshev+GMRES

% 'shf' shifts, 'wts' weights
% 'im'  Krylov subspace dimension
% 'tol' [tol of outer problems, tol of inner problem]
% 'dt' step size, 'num' num of fixed pnt iters, 'q' order
FCI = struct('np',np,'shf',[],'wts',[],'im',20,'tol',[2e-1,2e-1],...
             'nblock',nblock,'method',method, 'bet',[],'num',[],'q',[]);

dmax = MAT.rho(2);

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

if(nblock == 1)
    FCI.bet = [MAT.rho(1)/2-1, MAT.rho(1)/2, dmax];
else
    FCI.bet = [-1, dmax/2 + sqrt(dmax^2/4 + MAT.rho(1)), dmax];
end

% set up fixed point iterations
FCI.num = zeros(size(FCI.wts));
FCI.q   = zeros(size(FCI.wts));
FCI.d   = zeros(size(FCI.wts));
rates   = zeros(size(FCI.wts));
for p = 1:np
    [FCI.q(p), FCI.d(p), rates(p)] = exp_rate(FCI.shf(p), FCI.bet, 1, 7);
    FCI.num(p) = ceil( log(FCI.tol(1) * abs(FCI.wts(1)/FCI.wts(p)).^0.5 )/log(rates(p)));
end
FCI.num = max(FCI.num, 1);

rates
end