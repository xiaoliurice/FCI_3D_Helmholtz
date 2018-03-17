function [u,nmv] = fcisol(f,z0,MAT,FCI)
% one step of contour integration solution
% because of the last step of GMRES, outside this function an iterative refinement suffices

pfci=tic;

% outer problem
u = zeros(size(f));
nmv = 0;
nrm0 = norm(f);
for p = 1:FCI.np
   % apply the shift
   z  = z0 + FCI.shf(p);
   ss = sign(imag(z));
   tmp = sqrt(z + 1j*MAT.D(:));
   MAT.atmp = real(tmp).^2;
   MAT.btmp = (ss*imag(tmp)) ./ real(tmp);

   
   % matrix-free solution of shifted problems
   [v,nmv1] = cmplxsplit_exp(f,MAT,ss,FCI.num(p),FCI.dt(p),4);
   nmv = nmv + nmv1;
   fprintf('|%.2e',norm(f-helmop(v,z,MAT))/nrm0*abs(FCI.wts(p)/FCI.wts(1)));
   
   u = u + FCI.wts(p)*v;
end
v = helmop(u,z0,MAT);
c = (v'*f)/(v'*v); % step size
u = c*u;
v = f-c*v; % new residual
fprintf('|step %.2e,%.2e',real(c),imag(c));

% inner problem
% GMRES on current residual
[v,~,~,iter] = gmres(@(x)helmop(x,z0,MAT),v,FCI.im,FCI.tol(2),FCI.nim);
u = u + v;
nmv = nmv + (iter(1)-1)*FCI.im + iter(2) + 1;

t=toc(pfci); fprintf('|time %.2fs\n',t);
end