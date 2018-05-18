function [y] = helmop(x,MAT)
% Helmholtz operator
    n = prod(MAT.N);
    if(n==length(x))
        y = (stiffop(x,MAT) - (1+1j*MAT.D(:)).*MAT.M(:).*x)/(MAT.khmax^2);
    else
        y = [1j*x(n+1:end); (1j*stiffop(x(1:n),MAT) + (1+1j*MAT.D(:)).*MAT.M(:).*x(n+1:end)) / (-MAT.khmax^2)];
    end
end

function [y] = stiffop(x,MAT)
% matrix-vector product of the stiffness matrix
    if(MAT.sparse)
        y = MAT.S*x;
    else
        y = fftn(reshape(x,MAT.N));
        y(:) = MAT.S(:).*y(:);
        y = reshape(ifftn(y),[length(x),1]);
    end
end