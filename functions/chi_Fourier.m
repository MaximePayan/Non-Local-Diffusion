function chi = chi_Fourier(A,B,a,b,N,ipi)
% Computes the Fourier coefficients chi_0,chi_1,...,chi_N of the indicator 
% function of [a,b], defined on the interval [A,B] and extended as an 
% "even" periodic function. The normalization condition for the
% coefficients is such that
% chi(x) = chi_0 + 2 \sum_{n\geq 1} chi_n cos(n*pi*(x-A)/(B-A))
if nargin < 6
    ipi = pi;
end
n = (1:N-1)';
chi = [ (b-a)/(B-A);
        1./(n*ipi) .* ( sin(n*ipi*(b-A)/(B-A)) - sin(n*ipi*(a-A)/(B-A)) ) ];