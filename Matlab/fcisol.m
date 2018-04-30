function [u,nmv] = fcisol(f,MAT,FCI)
% one step of contour integration solution
% because of the last step of GMRES, outside this function an iterative refinement suffices

ftimer = tic;
% outer problem
u = zeros(size(f));
nmv = 0;
nrm0 = norm(f);
for p = 1:FCI.np
   z  = 1 + FCI.shf(p);
   
   % polynomial fixed point iteration of shifted problems
   %[v,nmv1] = exp_poly(f,MAT,z,FCI.num(p),FCI.dt(p),FCI.q);
   [v,nmv1] = cheby_poly(f,MAT,z,FCI.tol(2),FCI.num(p)*FCI.q);
   nmv = nmv + nmv1;
   fprintf('|%.2e',norm(f-helmop(v,z,MAT))/nrm0*abs(FCI.wts(p)/FCI.wts(1)));
   
   u = u + FCI.wts(p)*v;
end
v = helmop(u,1,MAT);
c = (v'*f)/(v'*v); % step size
u = c*u;
v = f-c*v; % new residual
fprintf('|outer step (%.2f,%.2f)',real(c),imag(c));

% inner problem
% GMRES on current residual
[v,~,relres,iter] = gmres(@(x)helmop(x,1,MAT),v,FCI.nim(2),FCI.tol(2),FCI.nim(1));
u = u + v;
nmv = nmv + (iter(1)-1)*FCI.nim(2) + iter(2) + 1;
fprintf('|inner res %.2e',relres);

t=toc(ftimer); fprintf('|time %.2fs\n',t);
end