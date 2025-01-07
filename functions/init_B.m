function [B,C,q] = init_B(which_B,ipi,N,Omega,arg1,arg2)
% B size N x N
Xmax = Omega(2);
Xmin = Omega(1);
l = (Xmax-Xmin);

if class(which_B) == "string" && (which_B == "Indicator" || "Heaviside" || "square")
    which_B = 1;
elseif class(which_B) == "string"
    which_B = 0;
end

if which_B
    if class(arg1)~= "cell" && class(arg2)~= "cell"
        Omega1 = arg1; %arg1 is interpreted as the indicator function on domain Omega1
        Omega2 = arg2;
        l1 = (sum(Omega1(:,2)) - sum(Omega1(:,1)));
        l2 = (sum(Omega2(:,2)) - sum(Omega2(:,1)));
        [nb_i1,~] = size(Omega1);
        [nb_i2,~] = size(Omega2);

        C = max([8*nb_i1*nb_i2*l/(ipi^2*l1),4*nb_i2/ipi,2*nb_i1*l2/(ipi*l1),l2/l]);
        
        % Omega1 and Omega2 are union of segments include in Omega
        chi1 = zeros(N,1);
        chi2 = zeros(N,1);
        for i1 = 1:nb_i1
            chi1 = chi1 + chi_Fourier(Xmin,Xmax,Omega1(i1,1),Omega1(i1,2),N); % Fourier coefficients of \chi_{Omega_1}
        end
        for i2 = 1:nb_i2
            chi2 = chi2 + chi_Fourier(Xmin,Xmax,Omega2(i2,1),Omega2(i2,2),N); % Fourier coefficients of \chi_{Omega_1}
        end
        chi2(2:end) = chi2(2:end)*2; % chi2'*u will give |domain|^{-1} * \int_{domain2} u 
        B = chi1/l1 * l*chi2';
        q = 1;
    else
        error("If indicator function case, specify domains Omega1 and Omega2")
    end
else
    if class(arg1)=="cell" && class(arg2)=="cell"
        g1 = arg1{1}; gp1 = arg1{2}; gpp1 = arg1{3};
        g2 = arg2{1}; gp2 = arg2{2}; gpp2 = arg2{3};
        %I any other case arg1 and arg2 are interpreted as cell
        %containing three handle functions e.g. {g,g',g"}
        chi1 = cos_Fourier(g1,N,Xmin,Xmax);
        chi2 = cos_Fourier(g2,N,Xmin,Xmax);
        X = linspace(Xmin,Xmax,200);
        C1 = l/ipi^2 * (abs(gp1(Xmax))+max(abs(gpp1(X))));
        C2 = l/ipi^2 * (abs(gp2(Xmax))+max(abs(gpp2(X))));
        chi2(2:end) = chi2(2:end)*2;
        B = l*chi1*chi2';
        C = l*C1*C2;
        q = 2; %by default, but it can be better for C1,C2 and q if a study is done on them
    else
        error("If no specification, give g1 and g2 two handle functions")
    end
end
end

