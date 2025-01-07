function r = norme_alpha(u,alpha,mat)
[p,q] = size(u);
k = [1,2.*(1:max(p,q)-1).^alpha];
if nargin == 2 || mat == false 
    r = k(1:p)*abs(u);
else
    r = k(1:p)*abs(u(:,1));
    for i=2:q
        r = max(r,k(1:p)*abs(u(:,i))./k(i));
    end
end
end