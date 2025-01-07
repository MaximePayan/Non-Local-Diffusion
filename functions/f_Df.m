function [f,Df] = f_Df(U,B,para)

if class(U) == "xAppErr"
    U = U.App;
end

N = length(U)/2;

domain = para.domain;
l = domain(2)-domain(1);
a = para.a;
b = para.b;
c = para.c;
d = para.d;
delta = para.delta;
theta = para.theta;

Lap = laplacien(N,l);
I = eye(N);


Df = [theta*diag(Lap) + a*I,b*I;
    c*I + delta*B, diag(Lap)+ d*I];

f = Df*U;

end
