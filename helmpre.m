function [y] = helmpre(x,N,z,kh)
% constant-coefficient Helmholtz solution
% dst(y) = dst(x)./(|xi|^2 - z*kh^2)
% it is a good preconditioner when |Im(z)| > |Re(z)|
    y = dst1fft( reshape(x,N) );
    xi = kron(ones(1,N(2)*N(3)),(pi/(N(1)+1) * (1:N(1))).^2) ...
        +kron( (pi/(N(3)+1) * (1:N(3))).^2, ones(1,N(1)*N(2)) )...
        +kron(ones(1,N(3)), kron( (pi/(N(2)+1) * (1:N(2))).^2, ones(1,N(1))));
    y(:) = y(:)./(xi(:)-z*kh^2);
    y = reshape(dst1fft(y),[length(x),1]);
end