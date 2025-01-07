function r = norme_nu(u,nu,mat)
k = ksi(length(u),nu);
if nargin == 2 || mat == false 
    r = k'*abs(u);
else
    r = k'*abs(u(:,1));
    for i=2:length(u(1,:))
        r = max(r,k'*abs(u(:,i))./k(i));
    end
end