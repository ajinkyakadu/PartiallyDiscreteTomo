function [L] = smoothReg(n,mode)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

switch mode
    case 1
        e1 = ones(n(2),1);
        e2 = ones(n(1),1);
        Lx = spdiags([e1 -2*e1 e1],-1:1,n(2),n(2));
        Lx(1,:) = 0; Lx(end,:) = 0;
        Lz = spdiags([e2 -2*e2 e2],-1:1,n(1),n(1));
        Lz(1,:) = 0; Lz(end,:) = 0;
        Ix = speye(n(2));
        Iz = speye(n(1));
        L = kron(Ix,Lz) + kron(Lx,Iz);
    case 2
        L = speye(prod(n));
end

end

