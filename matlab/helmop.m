function [y] = helmop(x,z,MAT)
% Helmholtz operator
% y = [DST*DD*DST - (z+1j*BD)*M] x
% z0 complex coefficient
    N = MAT.N;
    y = dst1fft(reshape(x, N));
    y(:) = MAT.S(:).*y(:);
    y = reshape(dst1fft(y), [length(x),1]);
    y = y - (z + 1j*MAT.D(:)).*(MAT.M(:).*x);
end