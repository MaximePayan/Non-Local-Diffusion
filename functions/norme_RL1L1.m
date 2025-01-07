function r = norme_RL1L1(U,alpha,mat)
[p,q] = size(U);
Np = (p-1)/2;
k = [1,2.*(1:max(Np,q)-1).^alpha];
if nargin == 2 || mat == false 
    r = abs(U(1)) + k(1:Np)*abs(U(2:Np+1)) + k(1:Np)*abs(U(Np+2:p));
else
    r = abs(U(1,1)) + k(1:Np)*abs(U(2:Np+1,1)) + k(1:Np)*abs(U(Np+2:p,1));
    %r2 = abs(U(1,q+1)) + k(1:Np)*abs(U(2:Np+1,q+1)) + k(1:Np)*abs(U(Np+2:p,q+1));
    for i=2:q
        r = max(r,(abs(U(1,i)) + k(1:Np)*abs(U(2:Np+1,i)) + k(1:Np)*abs(U(Np+2:p,i)))./k(i));
        %r2 = max(r2,(abs(U(1,q+i)) + k(1:Np)*abs(U(2:Np+1,q+i)) + k(1:Np)*abs(U(Np+2:p,q+i)))./k(i));
    end
    %r = max(r1,r2);
end
end