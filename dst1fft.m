function [x] = dst1fft(x)
% 3D DST-I from FFT
% its inverse is itself AA = I
N = size(x);

% first dimension
n =  N(1);
c = -0.5j/sqrt(2*n+2);
for i3 = 1:N(3)
    for i2 = 1:N(2)
        tmp = zeros(2*(n+1),1);
        tmp(2:n+1)   =  x(:,i2,i3);
        tmp(n+3:end) = -tmp(n+1:-1:2);
        tmp = fft(tmp);
        x(:,i2,i3) = (tmp(end:-1:n+3)-tmp(2:n+1))*c;
    end
end

% second dimension
n = N(2);
c = -0.5j/sqrt(2*n+2);
for i3 = 1:N(3)
    for i1 = 1:N(1)
        tmp = zeros(2*(n+1),1);
        tmp(2:n+1)   =  x(i1,:,i3);
        tmp(n+3:end) = -tmp(n+1:-1:2);
        tmp = fft(tmp);
        x(i1,:,i3) = (tmp(end:-1:n+3)-tmp(2:n+1))*c;
    end
end

% third dimension
n = N(3);
c = -0.5j/sqrt(2*n+2);
for i2 = 1:N(2)
    for i1 = 1:N(1)
        tmp = zeros(2*(n+1),1);
        tmp(2:n+1)   =  x(i1,i2,:);
        tmp(n+3:end) = -tmp(n+1:-1:2);
        tmp = fft(tmp);
        x(i1,i2,:) = (tmp(end:-1:n+3)-tmp(2:n+1))*c;
    end
end

end