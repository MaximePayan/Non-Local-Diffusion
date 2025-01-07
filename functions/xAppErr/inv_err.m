function inv_x = inv_err(x,nu)
    inv_x = xAppErr;
    if class(x) ~= "xAppErr"
        inv_x.App = x;
        x = inv_x;
    end
    one = eye(len(x),1);
    inv_x.App = convomat(x.App)\one;
    err = norm_err(convo(inv_x.App,x.App)-one,nu);
    ninv_x = norm_err(inv_x,nu);
    if err+ninv_x*x.Err > 1
        disp('Not invertible, check error')
    end
    inv_x.Err = ninv_x*(err+ninv_x*x.Err)./(1-err-ninv_x*x.Err);
end