function [A, z] = normalize(A, dim)

if(nargin < 2)
    z = sum(A(:));
    z(z==0) = 1;
    A = A./z;
else
    z = sum(A, dim);
    z(z==0) = 1;
    A = bsxfun(@rdivide, A, z);
end
end
