function [Radius,Center] = FarawayCircle(n,Srho_Ndiag,MajSum_Ndiag,C,f,params)
%%% Give the radii bounds above 
Radius = [];
Center = [];

a=params.a;
b=params.b;
c=params.c;
d=params.d;
delta = params.delta;
theta = params.theta;
domain = params.domain;
l = domain(2) - domain(1);

for i=1:length(n)
eps = mod(n(i)-1,2);
ip = (n(i)-1 - eps)/2 +1;
ipi = pi;
if eps == 1
    if class(c)=="intval"
            ip = intval(ip);
            ipi = intval('pi');
    end
    Radius =[Radius; C*delta*(Srho_Ndiag + MajSum_Ndiag).*f(ip)./ip + abs(c)] ;
    Center = [Center;-(ipi*ip/l)^2 + d];
else
    Radius = [Radius; abs(b)];
    Center = [Center;-theta*(ipi*ip/l)^2 + a];
end
end
