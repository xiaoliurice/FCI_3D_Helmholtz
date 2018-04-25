function [mdt, mrate] = exp_rate(rho,z,p)
   % search for the best step size and convergence rate
   % rho spectral radius of laplacian
   % z: complex factor in A-zM
   % p: order of approximation
   zrt = sqrt(z);
   rrt = abs(imag(zrt))/real(zrt); %|s|/r
   irt = sqrt(rho)/real(zrt);      % rho/r
   dt0 = rrt/(rrt*rrt+ irt*irt);   % a trial dt from min|(1-rrt*dt,irt*dt)|
   
   nlam = 20; % number of samples in imaginary parts
   mrate = 1;
   mdt = 0;
   for i = 1:10000
       lam = [complex(-rrt*ones(1,nlam), irt*(1:nlam)/nlam), 1j];
       dt = dt0*i;
       rate = convrate(dt, lam, p);
       if(rate < mrate)
          mdt = dt;
          mrate = rate;
       end
       if(rate > 1)
           return;
       end
   end
end

function [rate] = convrate(dt, z, p)
    c = ones(size(z));
    k = ones(size(z));
    for j = 1:p
        k = k.*z*(dt/j); %(zdt)^j/j!
        c = c + k;
    end
    
    rate = max(abs(c(1:end-1))/abs(c(end)));
end