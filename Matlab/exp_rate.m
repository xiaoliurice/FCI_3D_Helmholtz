function [z0, md, mrate] = exp_rate(z,rho,q,nblock)
   % search for the best step size and convergence rate
   % rho spectral radius of the Hermitian part and skew Hermitian part of S
   % z: complex B-zI
   % q: order of approximation
   
   r1 = rho(1);
   r2 = rho(2);
   zre = real(z);
   zim = imag(z);
   
   if(nblock == 1)
       b = [0,r1]-1-zre;
   else
       tmp = r2/2 + sqrt(r2*r2/4 + r1);
       b = [-tmp,tmp]-1-zre;
   end
   z0 = sum(b)/2 + 1j*zim; % z-z0 real
   db = (b(2) - b(1))/2;
   if(zim > 0 )
       b = b - 1j*zim;
   else
       b = b - 1j*(zim + r2);
   end
   
   ba = abs(b);
   io = imag(sum(ba)/sum(ba./b)+z); %imag(o): center-shift
   d0 = -0.25/(io-imag(z0));        %z0-i/d0 = o
   
   
   nlam = 20; % number of samples in imaginary parts
   mrate = 1;
   md = 0;
   lam = [complex(-zim*ones(1,nlam),db*(1:nlam)/nlam),...
          complex(-(zim+r2)*ones(1,nlam),db*(1:nlam)/nlam), -1j*(z-z0)];
   for i = 1:10000
       d = d0*i;
       rate = convrate(d, lam, q);
       if(rate < mrate)
          md = d;
          mrate = rate;
       end
       if(rate > 1)
           return;
       end
   end
   
end

function [rate] = convrate(d, z, q)
    c = ones(size(z));
    k = ones(size(z));
    for j = 1:q
        k = k.*z*(d/j);
        c = c + k;
    end
    
    rate = max(abs(c(1:end-1))/abs(c(end)));
end