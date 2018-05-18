function [d, rate] = ric_rate(z,rho,nblock)

r1 = rho(1);
r2 = rho(2);
zre = real(z);
zim = imag(z);

if(nblock == 1)
    b = [0,r1]-1-zre;
else
    tmp = r2/2 + sqrt(r2*r2/4 + r1);
    b = [-tmp,tmp]-1-zre;
end

if(zim > 0 )
    b = b - 1j*zim;
else
    b = b - 1j*(zim + r2);
end

ba = abs(b);
d  = sum(ba./b) / sum(ba);
rate = abs(b(2)-b(1)) / sum(ba);
end