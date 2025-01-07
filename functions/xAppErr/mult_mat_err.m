function M = mult_mat_err(u,m_output)
    if class(u)~="xAppErr"
        u = xAppErr(u);
    end
    M = convomat(u.App,m_output);
    M = xAppErr(M,u.Err);
end

