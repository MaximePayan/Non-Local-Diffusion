function r = norm_err(x,nu)
    if class(x)~="xAppErr"
        k = ksi(length(x),nu);
        r = k'*abs(x);
    else
        k= ksi(len(x),nu);
        r = k'*abs(x.App);
    end
end