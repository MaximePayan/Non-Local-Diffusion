function [A,W] = approx_inverse(A_bar,u,gv,gpv,Lap_inv,eps,N,D)

one = eye(2*N-1,1);
A21 = zeros(D,D);
if exist('intval','file')
    one = intval(one);
    A21 = intval(A21);
end

w11 = convomat(gv.App,2*N-1)\one; %size 2N-1
w12 = -convo(w11,convo(gpv.App,u,2*N-1))/eps;
w21 = zeros(2*N-1,1);
w22 = one/eps;

W = {w11, w12; w21, w22};

A11 = convomat(w11,D)*Lap_inv; %size D*D
A11(1:2*N-1,1:2*N-1) = A_bar(1:2*N-1,1:2*N-1);
A12 = convomat(w12,D)*Lap_inv;    
A12(1:2*N-1,1:2*N-1) = A_bar(1:2*N-1,2*N:4*N-2);

A21(1:2*N-1,1:2*N-1) = A_bar(2*N:4*N-2,1:2*N-1);
A22 = Lap_inv/eps;
A22(1:2*N-1,1:2*N-1) = A_bar(2*N:4*N-2,2*N:4*N-2);

A = {A11, A12; A21, A22};
end