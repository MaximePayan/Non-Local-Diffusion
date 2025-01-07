function m = mult_err(x,y,nu,trunc)
    m = xAppErr;
    if class(x) ~= "xAppErr"
        m.App = x;
        x = m;
    end
    if class(y) ~= "xAppErr"
        m.App = y;
        y = m;
    end
    m.App = convo(x.App,y.App,len(x)+len(y)-1);
    nx = norm_err(x.App,nu);
    ny = norm_err(y.App,nu);
    m.Err = nx*y.Err + ny*x.Err + x.Err*y.Err;
    if nargin == 4 && trunc<=len(m)
        e_trunc = norm_err(m.App - [m.App(1:trunc);zeros(len(m)-trunc,1)],nu);
        m.App = m.App(1:trunc);
        m.Err = m.Err+e_trunc;
    end
end