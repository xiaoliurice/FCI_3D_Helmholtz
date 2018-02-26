np   = 4;    % num of poles
rad1 = 3;    % in x
rad2 = 3;    % in y
theta = pi/np.*(1:2:np*2-1);  % angles in radi
shf = complex( rad1*cos(theta), rad2*sin(theta) ); % points on the ellipse
shf = shf - max(real(shf));                        % shift the center
wts = complex( rad2*cos(theta), rad1*sin(theta) )./shf/np;

for rnd = 1
    N = ones(1,3)*rnd*100; %grid
    kh0 = 2*pi/2.25;       %2pi/number of points per wavelength
    kh  = ones(N)*kh0;
    kh(ceil(N(1)/3):ceil(N(1)*2/3),ceil(N(2)/3):ceil(N(2)*2/3),ceil(N(3)/3):ceil(N(3)*2/3)) = kh0/2;
    
    % absorbing layers
    nab = 10;
    ab  = ones(N)*(kh0^2*1j);
    ab(nab+1:end-nab-1,nab+1:end-nab-1,nab+1:end-nab-1) = 0;
    ab = smooth3(ab,'gaussian',5);
    
    f = zeros(N);
    f(ceil(N(1)/2),ceil(N(2)/2),ceil(N(3)/2)) = 1;
    f = reshape(f,[prod(N),1]);
    nrm0 = norm(f(:));
    
    tol = 1e-6;
    u = zeros(size(f));
    for i = 1:60
        v = fcipre(f,N,kh,kh0,ab, shf,wts);
        u = u + v; 
        f = f - helmop(v,N,1,kh,ab);
        if (norm(f) < nrm0*tol)
            continue;
        end
    end
    
    u = reshape(u,N);
end