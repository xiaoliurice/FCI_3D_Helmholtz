function [mq, md, mrate] = exp_rate(z,bet,qmin,qmax)
   % search for the best step size and convergence rate
   % rho spectral radius of the Hermitian part and skew Hermitian part of S
   % z: complex B-zI
   % bet: spectrum info of B
   % qmin,qmax: min and max order of approximation
   
   zre = real(z);
   zim = imag(z);
   
   
   % optimal stat richardson
   b = [bet(1)-bet(2),bet(1)+bet(2)]-z;
   if(zim<0)
       b = b - 1j*bet(3);
   end
   ba = abs(b);
   md = -1/imag(sum(ba)/sum(ba./b));
   mrate = abs(b(2)-b(1)) / sum(ba);
   mq = 1;
   
   
   nlam = 20; % number of samples in imaginary parts
   lam = [complex(    -zim     *ones(1,nlam),bet(2)*(1:nlam)/nlam),...
          complex(-(zim+bet(3))*ones(1,nlam),bet(2)*(1:nlam)/nlam),-1j*(zre-bet(1))];
   
   d0 = md/10;
   for q = qmin:qmax
       for i = 1:10000
           d = d0*i;
           rate = convrate(d, lam, q)^(1/q);
           if(rate < mrate)
               mq = q;
               md = d;
               mrate = rate;
           end
           if(rate > 1)
               break;
           end
       end
   end
   mrate = mrate^mq;
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