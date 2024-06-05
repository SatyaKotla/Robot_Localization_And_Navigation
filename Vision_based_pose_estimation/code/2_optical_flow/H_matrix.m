function [H_i] = H_matrix(x,y,z)
     n = numel(x);
     Ap = [-1./z, zeros(n, 1), x./z;zeros(n, 1), -1./z, y./z];
     Bp = [x.*y, -(1+x.^2), y;(1+y.^2), -x.*y, -x];
     H_i = [Ap Bp];
end
     
