function fx=exp_err(alpha,x,K,nu,trunc)
%%% f(.) = exp(alpha . )
    fx = xAppErr;
    if class(x) ~= "xAppErr"
        fx.App = x;
        x = fx;
    end
    norm_x = norm_err(x,nu);
    abs_alpha = abs(alpha);
    if x.Err >0
        fx.Err = abs_alpha.*x.Err.*exp(abs_alpha*(norm_x+x.Err));
        x.Err = 0;
    end
    sum_err = exp(abs_alpha.*norm_x);
    N = len(x);
    xpoweri = alpha.*x;
    if nargin==5
        fx.App = eye(trunc,1);
        fact_i = 1;
        for i=1:K
            fact_i = i*fact_i;
            fx = fx + xpoweri./fact_i;
            sum_err = sum_err*abs_alpha.*norm_x/i;
            xpoweri = mult_err(alpha.*x,xpoweri,nu,trunc);
        end
    else
        fx.App = eye(K*(N-1)+1,1);
        fact_i = 1;
        for i=1:K
            fact_i = i*fact_i;
            fx = fx + xpoweri./fact_i;
            sum_err = sum_err*abs_alpha.*norm_x/i;
            xpoweri = mult_err(alpha.*x,xpoweri,nu);
        end
    end
    fx.Err = fx.Err + sum_err;
end