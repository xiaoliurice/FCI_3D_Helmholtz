function u = fcipre(f,N,kh,kh0,ab, shf,wts)
% Contour integration preconditioner
% shf vector of complex shifts
% wts vector of quadrature weights

restart = 20;
maxit   = 3;
tol = 1e-2;

tic;
% outer problem
u = zeros(size(f));
for p = 1:length(shf)
   z = 1 + shf(p);
   % matrix-free solution of shifted problems
   [v,~,relres,iter] = gmres(@(x)helmop(x,N,z,kh,ab),f,restart,tol,maxit,@(x)helmpre(x,N,z,kh0));
   fprintf('%dth pole: relative residual %.2e, iterations %d.\n',p,relres,(iter(1)-1)*restart + iter(2));
   u = u + wts(p)*v;
end

% inner problem unpreconditioned GMRES on current residual
v = f-helmop(u,N,1,kh,ab);
nrm = norm(v);
[v,~,relres] = gmres(@(x)helmop(x,N,1,kh,ab),v,restart,tol,maxit);
u = u + v;

toc;
fprintf('res %.2e, inner relres %.2e\n', nrm*relres,relres); % report current residual
end