function [Upp,E] = func_Newton(U,F_DF,it_max,tol_min,tol_max)
    if nargin<3
        tol_min=10^-12;
        tol_max=10^6;
        it_max=20;
    end
    itt = 0;
    Upp=U;
    [FU,DFU] = F_DF(Upp);
    e=norm(FU,1);
    E=e;
    while e>tol_min && e<tol_max && itt<it_max
        Upp = Upp - DFU\FU;
        [FU,DFU]=F_DF(Upp);
        e = norm(FU,1);
        E=[E,e];
        itt = itt + 1;
    end
end

