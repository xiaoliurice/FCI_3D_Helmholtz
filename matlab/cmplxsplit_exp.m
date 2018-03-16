function [u, nmv] = cmplxsplit_exp(f,MAT,ss,num,dt,p)
% solution of shifted Helmholtz via p-th order taylor expansion of exp func

u  = zeros(size(f));
uu = zeros(size(f));

% compensation factors
w   = -1j*ss;
wdt = w*dt;
wf = zeros(1,p);
t = 1/w;
for j=1:p
    t = t*wdt/j;
    wf(j) = t; %wdt^j/j!/w
end
wa = 1+w*sum(wf);

for i = 1:num
    k  = u;
    kk = uu;
    for j=1:p
        [k,kk] = sysop(k*(dt/j),kk*(dt/j),f,wf(j),MAT);
        u  = u  + k;
        uu = uu + kk;
    end
    u  =  u/wa;
    uu = uu/wa;
end
nmv = num*p;
end