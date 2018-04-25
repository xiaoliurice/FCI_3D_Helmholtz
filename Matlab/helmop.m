function [y] = helmop(x,z,MAT)
% Helmholtz operator
% y = [stiffness - (z+1j*MAT.D)*MAT.M]*x
% z complex shift
    y = stiffop(x,MAT) - (z+1j*MAT.D(:)).*(MAT.M(:).*x);
end