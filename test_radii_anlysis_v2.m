function [RadiiBounds,prec,M4,faraway_max]=test_radii_anlysis_v2(params,ipi,parameter_func,parameter_B,prec)

if nargin == 4 && class(ipi) == "intval"
    error("need argument prec, for precalculation without intval")
end

if nargin == 5 && class(ipi) ~= "intval"
    warning("calculations are made without intlab")
end

delta = params.delta;
theta = params.theta;
a = params.a;
b = params.b;
c = params.c;
d = params.d;
Omega = params.domain;

Xmax = Omega(2);
Xmin = Omega(1);
l = (Xmax-Xmin);

f = parameter_func{1}; 
p = parameter_func{2};
which_f = parameter_func{3}; %exp(px)/x, logx^p, x^p

B = parameter_B{1};
C = parameter_B{2}; %B s.t. |B_ij| <= C/(i^q*j^q)
q = parameter_B{3};

Ndiag = size(B,1);
% bound1 = Ndiag; to use when f is not a power
% bound2 = +Inf;
NK = [1, (1:Ndiag-1)];
RadiiBounds = zeros(2*Ndiag,1);

if class(ipi)=="intval"
    NK = intval(NK);
%    bound1 = midrad(mid(bound1),0);
%    bound2 = midrad(mid(bound2),0);
%    B = intval(B);
    C = intval(C);
    q = intval(q);
    RadiiBounds = intval(RadiiBounds);
    p = intval(p);
end

lambda = -(ipi/l*(0:Ndiag-1)).^2; % eigenvalues of the laplacian

%% The matrix of the linearized operator, but in "our" basis
M1 = [diag(theta*lambda+a), b*eye(Ndiag,Ndiag);
          c*eye(Ndiag,Ndiag)+delta*B, diag(lambda+d)];
%% Change of basis toward "your" basis
S = zeros(2*Ndiag,2*Ndiag,'logical');
for n = 1:Ndiag
    S(1+2*(n-1),n) = 1; 
    S(2+2*(n-1),Ndiag+n) = 1;
end
Sinv = S';

M2 = S*M1*Sinv;

%% Change of basis to make each row summable
Q = diag(repelem([1, 1:Ndiag-1],2).^p);
M3 = (Q*M2)/Q;
%% Diagonalization of the first few modes
if class(ipi)=="intval"
    Pprec = prec{1};
    Pprec_inv = prec{2};
    Dprec = prec{3};
else
    [Pprec,Dprec] = my_eig(M3);
    Pprec_inv = inv(Pprec);
    prec = {Pprec,Pprec_inv,Dprec};
end

M4 = Pprec\M3*Pprec;

%% Radii Bounds

switch which_f
    case "exp"
        MajSum_Ndiag  = exp(-p*Ndiag)/((1 - exp(-p))*Ndiag^(q-1));
    case "log"
        MajSum_Ndiag = 1./((p-1).*log(Ndiag).^(p-1).*Ndiag^(q-1));
    case "pow"
        MajSum_Ndiag = 1./((p+q-1).*Ndiag.^(p+q-1));
    otherwise
        if length(parameter_func) > 3
            MajSum_Ndiag = parameter_func{4};
        else
            error("You need to specify a bound on the sum_{j=Ndiag}^{+\infty} 1/(f(j)j)")
        end
end

rho_0 = zeros(Ndiag,1);
rho_1 = zeros(Ndiag,1);
if class(ipi) == "intval"
    rho_0 = intval(rho_0);
    rho_1 = intval(rho_1);
end
for i=1:Ndiag
    rho_i0 = f(NK).*abs(Pprec_inv(2*i-1,2.*(1:Ndiag)))./NK.^q;
    rho_i1 = f(NK).*abs(Pprec_inv(2*i,2.*(1:Ndiag)))./NK.^q;
    rho_0(i) = sum(rho_i0);
    rho_1(i) = sum(rho_i1); 
    RadiiBounds(2*i-1) = C*delta*rho_0(i)*MajSum_Ndiag + sum(abs(M4(2*i-1,[1:2*i-2 2*i:2*Ndiag])));
    RadiiBounds(2*i) = C*delta*rho_1(i)*MajSum_Ndiag + sum(abs(M4(2*i,[1:2*i-1 2*i+1:2*Ndiag])));
end

Srho_Ndiag = 0;
for k=1:Ndiag
    for j=1:Ndiag
        kp = max([1,k-1]);
        if class(c)=="intval"
            kp = intval(kp);
        end
        Srho_Ndiag = Srho_Ndiag+abs(Pprec(2*k-1,2*j-1) + Pprec(2*k-1,2*j))./(f(kp).*kp.^q);
    end
end

prec = {Pprec,Pprec_inv,Dprec,rho_0,rho_1,C*delta*Srho_Ndiag,C*delta*MajSum_Ndiag};


switch which_f
    case "pow"
        faraway_max = +Inf;
        b0 = abs(b) -theta*(ipi*Ndiag/l)^2 + a;
        b1 = C*delta*(Srho_Ndiag + MajSum_Ndiag).*f(Ndiag)./Ndiag.^q + abs(c) -(ipi*Ndiag/l)^2 + d;
        if (p > 0 && q > 0 && p - q <= 0) || (p - q == 2 && C*delta*Srho_Ndiag - (ipi/l)^2 <= 0)
            faraway_max = max([b0,b1]);
        elseif p - q > 0 && p - q < 2
            x_m = (2*ipi^2/(C*delta*l^2*Srho_Ndiag*(p-q))).^(1/(p-q-2));
            if class(x_m) == "intval"
                i_m = floor(x_m.inf):ceil(x_m.sup);
                i_m = max(Ndiag,i_m);
                bm  = C*delta*(Srho_Ndiag + MajSum_Ndiag).*f(i_m)./i_m.^q + abs(c) -(ipi.*i_m./l).^2 + d;
                faraway_max = max([bm, b0]);
            else
                i_m = floor(x_m):ceil(x_m);
                i_m = max(Ndiag,i_m);
                bm  = C*delta*(Srho_Ndiag + MajSum_Ndiag).*f(i_m)./i_m.^q + abs(c) -(ipi.*i_m./l).^2 + d;
                faraway_max = max([bm, b0]);
            end
        end
    case "exp"
        faraway_max = +Inf; %indeed this is not conclusive ...
    case "log"
        faraway_max = +Inf;
        x_m = func_Newton(Ndiag,@f_fp);
        b0 = abs(b) -theta*(ipi*Ndiag/l)^2 + a;
        if class(x_m) == "intval"
            i_m = floor(x_m.inf):ceil(x_m.sup);
            i_m = max(Ndiag,i_m);
            bm  = C*delta*(Srho_Ndiag + MajSum_Ndiag).*f(i_m)./i_m.^q + abs(c) -(ipi.*i_m./l).^2 + d;
            faraway_max = max([bm, b0]);
        else
            i_m = floor(x_m):ceil(x_m);
            i_m = max(Ndiag,i_m);
            bm  = C*delta*(Srho_Ndiag + MajSum_Ndiag).*f(i_m)./i_m.^q + abs(c) -(ipi.*i_m./l).^2 + d;
            faraway_max = max([bm, b0]);
        end
    otherwise
        faraway_max = nan;
        warning("No bound on faraway circles is computed since f is not in the study")
end
% for f = log^p
function [fx,fpx] = f_fp(x)
    fx = C*delta*(Srho_Ndiag + MajSum_Ndiag).*log(x).^p./x.^q + abs(c) -(ipi.*x./l).^2 + d;
    fpx = C*delta*(Srho_Ndiag + MajSum_Ndiag).*log(x).^(p-1).*(p-q.*log(x))./x.^(q+1)-2.*(ipi./l).^2.*x;
end
end
