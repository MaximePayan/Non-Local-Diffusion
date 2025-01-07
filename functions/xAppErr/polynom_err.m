function Px = polynom_err(x,P,nu,trunc)
Px = xAppErr;
if class(x) ~= "xAppErr"
    Px.App = x;
    x = Px;
end
if class(P)~="Polynomial"
    P = Polynomial(P);
end
degP =deg(P);
N = len(x);
% if degP <= 0
%     Px.App = scalar_eval(0,P).*eye(N,1);
%     Px.Err = 0;
% else
    if x.Err > 0
        Px.Err = x.Err.*scalar_eval(norm_err(x,nu)+x.Err,abs(derivative(P)));
        x.Err = 0;
    end
    xpower = 1;
    if nargin==4
        Px.App = P.Value(1);
        for i=1:degP
            xpower = mult_err(x,xpower,nu,trunc);
            Px = Px + P.Value(i+1).*xpower;
        end
    else
        Px.App = P.Value(1);
        for i=1:degP
            xpower = mult_err(x,xpower,nu);
            Px = Px + P.Value(i+1).*xpower;
        end
    end
% end
end