function [u,nmv] = fcisol2(f,MAT,FCI)
% one step of contour integration solution
% because of the last step of GMRES, outside this function an iterative refinement suffices

ftimer = tic;
% outer problem
u = zeros(size(f));
nmv = 0;
nrm0 = norm(f);
for p = 1:FCI.np
    z = FCI.shf(p);
    tol = FCI.tol(1) * sqrt(abs(FCI.wts(1)/FCI.wts(p)));
    
   % polynomial fixed point iteration of shifted problems
    [v,nmv1, relres] =   exp_poly(f,z,MAT,FCI.nblock, tol, FCI.num(p), FCI.z0(p),FCI.d(p),FCI.q);
  %  [v,nmv1, relres] = cheby_poly(f,z,MAT,FCI.nblock, tol, FCI.num(p)*FCI.q);
    fprintf('|%.1e,%d',relres*abs(FCI.wts(p))/abs(FCI.wts(1)),nmv1);
    nmv = nmv + nmv1;
   
  %  w = f - helmshift(v,z,MAT);
  %  nrm = norm(w);
  %  if(nrm > nrm0*tol)
  %     [w,~,relres,iter] = gmres(@(x)helmshift(x,z,MAT), w, FCI.im, nrm0*tol/nrm, FCI.num(p)*FCI.q);
  %     v = v + w;
  %     nmv1 = (iter(1)-1)*FCI.im + iter(2);
  %     fprintf(';%.1e,%d', nrm*relres/nrm0 * abs(FCI.wts(p))/abs(FCI.wts(1)), nmv1);
  %     nmv = nmv + nmv1;
  %  end
    
    u = u + FCI.wts(p)*v;
end
v = helmop(u,MAT);
c = (v'*f)/(v'*v); % step size
u = c*u;
v = f-c*v; % new residual
fprintf('|outer step (%.2f,%.2f)',real(c),imag(c));

% inner problem
% GMRES on current residual
[v,~,relres,iter] = gmres(@(x)helmop(x,MAT),v,FCI.im,FCI.tol(2),ceil(nmv1/FCI.im));
u = u + v;
nmv1 = (iter(1)-1)*FCI.im + iter(2) + 1;
nmv = nmv + nmv1;
fprintf('|inner %.2e, %d',relres, nmv1);

t=toc(ftimer); fprintf('|time %.2fs\n',t);
end

function [y] = helmshift(x,z,MAT)
    y = helmop(x,MAT) - z*x;
end

