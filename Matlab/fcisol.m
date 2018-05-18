function [u,nmv] = fcisol(f,MAT,FCI)
% one step of contour integration solution
% because of the last step of GMRES, outside this function an iterative refinement suffices

ftimer = tic;
% outer problem
u = zeros(size(f));
nmv = 0;
for p = 1:FCI.np
   
   % polynomial fixed point iteration of shifted problems=
   [v,nmv1, relres] = exp_poly(f,FCI.shf(p),MAT,FCI.nblock, FCI.tol(1), FCI.num(p), FCI.z0(p),FCI.d(p),FCI.q);
   fprintf('|%.1e',relres);
   nmv = nmv + nmv1;
   
   u = u + FCI.wts(p)*v;
end
v = helmop(u,MAT);
c = (v'*f)/(v'*v); % step size
u = c*u;
v = f-c*v; % new residual
fprintf('|outer step (%.2f,%.2f)',real(c),imag(c));

% inner problem
% GMRES on current residual
[v,~,relres,iter] = gmres(@(x)helmop(x,MAT),v,FCI.nim(2),FCI.tol(2),FCI.nim(1));
u = u + v;
nmv = nmv + (iter(1)-1)*FCI.nim(2) + iter(2) + 1;
fprintf('|inner res %.2e',relres);

t=toc(ftimer); fprintf('|time %.2fs\n',t);
end