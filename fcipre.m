function u = fcipre(f,MAT,FCI)
% Contour integration preconditioner

% outer problem
u = zeros(size(f));

tic;
for p = 1:FCI.np
   MAT.z = 1 + FCI.shf(p);
   % matrix-free solution of shifted problems
   [v,~,rres,iter] = gmres(@(x)helmop(x,MAT),f,FCI.im,FCI.tol,FCI.nim,@(x)helmpre(x,MAT));
   fprintf('|%d,%.2e',(iter(1)-1)*FCI.im + iter(2),rres);
   u = u + FCI.wts(p)*v;
end
fprintf('|\n');

% inner problem unpreconditioned GMRES on current residual
v = f-helmop(u,MAT);
nrm = norm(v);
MAT.z = 1;
[v,~,rres] = gmres(@(x)helmop(x,MAT),v,FCI.im,FCI.tol,FCI.nim);
u = u + v;
toc;

% report current residual
fprintf('\n res %.2e, inner rres %.2e\n', nrm*rres,rres); 
end