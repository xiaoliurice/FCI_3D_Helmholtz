function [u, nmv] = exp_poly(f,MAT,z,num,dt,q)
% solution of shifted Helmholtz via p-th order taylor expansion of exp func

ss = sign(imag(z));
tmp = sqrt(z + 1j*MAT.D(:));
atmp = real(tmp).^2;
btmp = ss*imag(tmp) ./ real(tmp); 

u  = zeros(size(f));
uu = zeros(size(f));

% compensation factors
w   = -1j*ss;
wdt = w*dt;
wf = zeros(1,q);
t = 1/w;
for j=1:q
    t = t*wdt/j;
    wf(j) = t; %wdt^j/j!/w
end
wa = 1+w*sum(wf);

for i = 1:num
    k  = u;
    kk = uu;
    for j=1:q
        [k,kk] = sysop(k*(dt/j),kk*(dt/j),f,wf(j),atmp,btmp,MAT);
        u  = u  + k;
        uu = uu + kk;
    end
    u  =  u/wa;
    uu = uu/wa;
end
nmv = num*q;
end

function [y1,y2] = sysop(x1,x2,f,w,atmp,btmp, MAT)
% 2x2 system operator
%  / -b     1 \ /x1\   /   0   \
%  \-M\S/a  -b/ \x2/ + \w*M\f/a/

y1 = x2 - btmp.*x1;
y2 = (w*f - stiffop(x1,MAT))./MAT.M(:)./atmp - btmp.*x2;
end