function xpowern = power_err(x,n,nu,trunc)
    xpowern = xAppErr;
    if class(x) ~= "xAppErr"
        xpowern.App = x;
        x = xpowern;
    end
    xpowern.App = x.App;
    xpowern.Err = x.Err;
    if nargin == 4
        if len(x) < trunc
            xpowern.App = zeros(trunc,1);
        end
        for i=1:len(x)
            xpowern.App(i) = x.App(i);
        end
        for i=2:n
            xpowern = mult_err(xpowern,x,nu,trunc);
        end
    else
        for i=2:n
            xpowern = mult_err(xpowern,x,nu);
        end
    end

end