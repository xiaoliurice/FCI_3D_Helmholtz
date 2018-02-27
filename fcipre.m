function [u,nmv] = fcipre(f,MAT,FCI)
% Contour integration preconditioner

% outer problem
u = zeros(size(f));
nmv = 0;
tic;
for p = 1:FCI.np
   MAT.z = 1 + FCI.shf(p);
   % matrix-free solution of shifted problems
   [v,~,rres,iter] = gmres(@(x)helmop(x,MAT),f,FCI.im,FCI.tol,FCI.nim,@(x)helmpre(x,MAT));
   fprintf('|%d,%.1e',(iter(1)-1)*FCI.im + iter(2),rres);
   nmv = nmv + 2*((iter(1)-1)*FCI.im + iter(2));
   u = u + FCI.wts(p)*v;
end

% inner problem unpreconditioned GMRES on current residual
MAT.z = 1;
v = f-helmop(u,MAT);
[v,~,~,iter] = gmres(@(x)helmop(x,MAT),v,FCI.im,FCI.tol,FCI.nim);
u = u + v;
nmv = nmv + (iter(1)-1)*FCI.im + iter(2) + 1;
t=toc;
fprintf('|%.2fs\n',t);
end