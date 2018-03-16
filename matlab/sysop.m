function [y1,y2] = sysop(x1,x2,f,w,MAT)
% 2x2 system operator
%  / -b     1 \ /x1\   /   0 \
%  \-M\S/a  -b/ \x2/ + \M\f/a/

y1 = x2 - MAT.btmp.*x1;
y2 = dst1fft(reshape(x1, MAT.N));
y2(:) = MAT.S(:).*y2(:);
y2 = reshape(dst1fft(y2), [length(x1),1]);
y2 = (w*f - y2)./MAT.M(:)./MAT.atmp -MAT.btmp.*x2;
end