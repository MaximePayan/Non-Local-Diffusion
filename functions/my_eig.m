function [V,D] = my_eig(A)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[V,D] = eig(A);

[d,ind] = sort(diag(D),'ComparisonMethod','real');

V = V(:,flip(ind));
D = diag(flip(d));

end

