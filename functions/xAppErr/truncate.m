function x_trunc = truncate(x,nu,trunc)
    x_trunc = xAppErr;
    if class(x)~="xAppErr"
        x_trunc.App = x;
        x = x_trunc;
    end
    x_trunc.Err = x.Err + norm_err([zeros(trunc,1);x.App(trunc+1:end)],nu);
    x_trunc.App = x.App(1:trunc);
end