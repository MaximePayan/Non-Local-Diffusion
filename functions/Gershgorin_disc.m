function [Center,Radius] = Gershgorin_disc(M,on_columns)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin<2
    on_columns = false;
end
Center = diag(M);
N = length(Center);
Radius = zeros(N,1);
if class(M) == "intval"
    Radius = intval(Radius);
end
if on_columns
    transposeM = M';
    for i=1:N
        Radius(i) = sum(abs(transposeM(i,[1:i-1,i+1:end])));
    end
else
    for i=1:N
        Radius(i) = sum(abs(M(i,[1:i-1,i+1:end])));
    end
end
end

