function [y] = helmop(x,N,z,kh,ab)
% Helmholtz operator
% y = -Laplacian(x) - (z*kh^2 + ab.)*x
% Laplacian via DST (implying zero boundary)
% z complex coefficient (1+shift)
% kh variable-coefficient dimensionless wavenumber  2pi/ppw in [0 -- pi (Nyquist rate)]
% ab variable-coefficient absorption in the form of 1j*kh^2*const
    y = dst1fft(reshape(x, N));
    xi = kron(ones(1,N(2)*N(3)),(pi/(N(1)+1) * (1:N(1))).^2) ...
        +kron( (pi/(N(3)+1) * (1:N(3))).^2, ones(1,N(1)*N(2)) )...
        +kron(ones(1,N(3)), kron( (pi/(N(2)+1) * (1:N(2))).^2, ones(1,N(1))));
    y(:) = y(:).*xi(:);
    y = reshape(dst1fft(y), [length(x),1]);
    y = y - (ab(:)+z*kh(:).^2).*x;
end