function [y] = helmop(x,MAT)
% Helmholtz operator
    n = prod(MAT.N);
    if(n==length(x))
        y = stiffop(x,MAT) - (MAT.M(:)+1j*MAT.D(:)).*x;
    else
        % double sized one
        y = [1j*x(n+1:end) - x(1:n); -1j*(stiffop(x(1:n),MAT)+(1-MAT.M(:)).*x(1:n)) - (1+1j*MAT.D(:)).*x(n+1:end)];
    end
end