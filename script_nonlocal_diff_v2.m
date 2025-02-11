function [A,E_maj,chi_q,bounds,prec] = script_nonlocal_diff_v2(parameter_B,params,D,alpha,prec)

if nargin == 4 && class(params.delta) == "intval"
    error("need argument prec, for precalculation without intval")
end

if nargin == 5 && class(params.delta) ~= "intval"
    warning("calculations are made without intlab")
end
B = parameter_B{1};
C = parameter_B{2};
q = parameter_B{3};

Ndiag = size(B,1);

if class(params.delta) == "intval"
    ipi = intval('pi');
    Ubartilde = prec{1};
    Utilde = prec{2};
    F = @(Utilde) functional(Utilde,Ubartilde,@(U) f_Df(U,B,params));
    [FU,DFU] = F(Utilde);
    A_app = prec{3};
    A = intval(zeros(2*(Ndiag+D)+1,2*(Ndiag+D)+1));
    I = intval(eye(2*Ndiag+1));
else
    ipi = pi;
    if nargin == 5 && class(prec) == "double"
        Ubartilde = prec;
    else
        [~,Df] = f_Df(zeros(2*Ndiag,1),B,params);
        [eigVec,eigVal] = my_eig(Df);
        Ubartilde = [eigVal(1,1); eigVec(:,1)];
    end
    F = @(Utilde) functional(Utilde,Ubartilde,@(U) f_Df(U,B,params));
    [Utilde,~] = func_Newton(Ubartilde,F);
    [FU,DFU] = F(Utilde);
    A_app = inv(DFU);
    prec = {Ubartilde,Utilde,A_app};
    A = zeros(2*(Ndiag+D)+1,2*(Ndiag+D)+1);
    I = eye(2*Ndiag+1);
end

l = params.domain(2)-params.domain(1);



%%
A(1,1) = A_app(1,1);
A(1,2:Ndiag+1) = A_app(1,2:Ndiag+1);
A(2:Ndiag+1,2:Ndiag+1) = A_app(2:Ndiag+1,2:Ndiag+1);
A(1,Ndiag+D+2:2*Ndiag+D+1) = A_app(1,Ndiag+2:end);
A(2:Ndiag+1,Ndiag+D+2:2*Ndiag+D+1) = A_app(2:Ndiag+1,Ndiag+2:end);
A(2:Ndiag+1,1) = A_app(2:Ndiag+1,1);
A(Ndiag+D+2:2*Ndiag+D+1,1) = A_app(Ndiag+2:end,1);
A(Ndiag+D+2:2*Ndiag+D+1,2:Ndiag+1) = A_app(Ndiag+2:end,2:Ndiag+1);
A(Ndiag+D+2:2*Ndiag+D+1,Ndiag+D+2:2*Ndiag+D+1) = A_app(Ndiag+2:end,Ndiag+2:end);
A(Ndiag+2:Ndiag+D+1,Ndiag+2:Ndiag+D+1) = diag((params.a-Utilde(1)-params.theta.*((Ndiag+1:Ndiag+D).*ipi/l).^2).^(-1));
A(2*Ndiag+D+2:2*(Ndiag+D)+1,2*Ndiag+D+2:2*(Ndiag+D)+1) = diag((params.d-Utilde(1)-((Ndiag+1:Ndiag+D).*ipi/l).^2).^(-1));
%% Y
E_maj = l^2/((1+q-alpha)*ipi^2)*(Ndiag-1-l/ipi*sqrt(max(0,params.d-Utilde(1))))^(alpha-(1+q));
chi_q = [1,1:Ndiag-1]'.^(-q);
R1_N = 1/abs(-params.theta*(l*Ndiag/ipi)^2 + (params.a-Utilde(1)));
R2_N = 1/abs(-(l*Ndiag/ipi)^2 + (params.d-Utilde(1)));

Y = norme_RL1L1(A([1:Ndiag+1,Ndiag+D+2:2*Ndiag+D+1],[1:Ndiag+1,Ndiag+D+2:2*Ndiag+D+1])*FU,alpha) + ...
    2*C*abs(params.delta)*abs(Utilde(2:Ndiag+1))'*chi_q*E_maj;
%% Z1
Mat = I - A([1:Ndiag+1,Ndiag+D+2:2*Ndiag+D+1],[1:Ndiag+1,Ndiag+D+2:2*Ndiag+D+1])*DFU;
Z1a = norme_RL1L1(Mat,alpha,true);
Z1b = C*abs(params.delta)/(2*Ndiag^(q+alpha)) * ((abs(A(1,Ndiag+D+2:2*Ndiag+D+1))*chi_q) + norme_alpha(abs(A(2:Ndiag+1,Ndiag+D+2:2*Ndiag+D+1))*chi_q,alpha) + norme_alpha(abs(A(Ndiag+D+2:2*Ndiag+D+1,Ndiag+D+2:2*Ndiag+D+1))*chi_q,alpha)) ...
    + abs(params.c)*R2_N+2*C*abs(params.delta)*E_maj;
Z1c = abs(params.b)*R1_N;
Z1 = max([Z1a, Z1b, Z1c]);
%% Z2
% Z2a = norme_RL1L1(A([1:Ndiag+1,Ndiag+D+2:2*Ndiag+D+1],1),alpha);
% Z2b = norme_RL1L1(A([1:Ndiag+1,Ndiag+D+2:2*Ndiag+D+1],2:Ndiag+1),alpha,true);
% Z2b = max([Z2b R1_N]);
% Z2c = norme_RL1L1(A([1:Ndiag+1,Ndiag+D+2:2*Ndiag+D+1],Ndiag+D+2:2*Ndiag+D+1),alpha,true);
% Z2c = max([Z2c R2_N]);
% Z2 = max([Z2a Z2b Z2c]);
Z2 = max([norme_RL1L1(A_app,alpha,true), R1_N, R2_N]);
bounds = [Y,Z1,Z2];
%% 
disc = (1-Z1)^2 - 4*Z2*Y;
if disc > 0
    r_min = ((1-Z1) - sqrt(disc))/(2*Z2); r_max = ((1-Z1) + sqrt(disc))/(2*Z2);
    if class(ipi) == "intval"
        r_min = r_min.sup;
        r_max = r_max.inf;
    end
    if r_max>0
        disp("There existe a unique solution in a r-neighbourhood where r belongs to ["+num2str(max(0,r_min))+","+num2str(r_max)+"].")
    else
        disp("We cannot conclude on the existence"); %Z1 must be smaller than 1
    end
else
    r_min = nan; r_max = nan;
    disp("We cannot conclude on the existence") %Z2 too large
end
bounds = [bounds, r_min, r_max];

end