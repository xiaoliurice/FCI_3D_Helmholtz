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