function [u,nmv] = fcisol(f,MAT,FCI)
% one step of contour integration solution
% because of the last step of GMRES, outside this function an iterative refinement suffices

pfci=tic;

% outer problem
u = zeros(size(f));
nmv = 0;
for p = 1:FCI.np
   MAT.zf = 1 + FCI.shf(p); % shift of forward operator
   MAT.zs = MAT.zf;         % shift for the solution operator
   % matrix-free solution of shifted problems
   [v,~,rres,iter] = gmres(@(x)helmop(x,MAT),f,FCI.im,FCI.tol,FCI.nim,@(x)helmpre(x,MAT));
   fprintf('|%d,%.1e',(iter(1)-1)*FCI.im + iter(2),rres);
   nmv = nmv + 2*((iter(1)-1)*FCI.im + iter(2));
   u = u + FCI.wts(p)*v;
end
MAT.zf = 1;
v = helmop(u,MAT);
c = real(v'*f) / (v'*v); % step size
u = c*u;
v = f-c*v; % new residual
fprintf('|%.2f',c);

% inner problem
% GMRES on current residual
MAT.zs = 1 + FCI.shf(1);
[v,~,~,iter] = gmres(@(x)helmop(x,MAT),v,FCI.im,FCI.tol,FCI.nim,@(x)helmpre(x,MAT));
u = u + v;
nmv = nmv + (iter(1)-1)*FCI.im + iter(2) + 1;

t=toc(pfci); fprintf('|%.2fs\n',t);
end