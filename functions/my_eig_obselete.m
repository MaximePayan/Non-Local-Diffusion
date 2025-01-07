function [V,D] = my_eig_obselete(A,vz,nb_eig,err,limit,eps)
%Compute the eigen elements of A by the power methode
%   A is a square matrix and U an unitary matrix for initialisation
% D contains the eigvalues, V(i) is the eigvector of Ai associated to D(i)
% (execpt i=0, V(i) is not a eigvector of A ...)
N = length(A);
if nargin == 1
    vz = eye(N,1);
end
if nargin < 3
    err = 1e-10;
    limit = 1e4;
    nb_eig = 1;
    eps = 0.1;
end

vz = loop_my_eig(A,vz,err,limit);
z = (A*vz)./vz;

if abs(z - sum(z)/N) > err %If no convergence then shift
    B = A + eps.*(1+1i)*eye(N);
    vz = loop_my_eig(B,vz,err,limit);
    z = (B*vz)./vz - eps.*(1+1i);
end
V = vz;
D = z;
% For the other eigenelements
Ai = A;
for i = 2:nb_eig
    Ai = Ai - z.*vz*vz';
    if nargin == 1
        vz = zeros(N,1);
        vz(i) = 1;
    end

    vz = loop_my_eig(Ai,vz,err,limit);
    z = (Ai*vz)./vz;

    if abs(z - sum(z)/N) > err %If no convergence then shift
        B = Ai + eps.*(1+1i)*eye(N);
        vz = loop_my_eig(B,vz,err,limit);
        z = (Ai*vz)./vz;
    end
    V = [V,vz];
    D = [D,z];
end



end

function v = loop_my_eig(B,v0,err,limit)
    K = length(B);
    vp = zeros(K,1);
    v = v0;
    j = 0;
    while norm(v - vp) > err && j<limit
        vp = v;
        v = B\vp;
        v = v/norm(v);
        j=j+1;
    end
end
