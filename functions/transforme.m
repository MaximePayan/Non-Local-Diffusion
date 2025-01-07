function [y,y_p] = transforme(w,a,b,g,g_prime)
% w is a vector contains fourier's coeff, return fourier's coeff of
% gamma(w_0 + 2*\sum w_k*cos(k*pi/l . ))
N =length(w);
l=b-a;

y = zeros(N,1);
y_p = zeros(N,1);
func_w = @(x) g(Fourier_cos(w,x,l));
func_p_w = @(x) g_prime(Fourier_cos(w,x,l));
y(1) = integral(func_w,a,b)/l;
y_p(1) = integral(func_p_w,a,b)/l;
for n=2:N
    fw = @(x) func_w(x).*cos((n-1)*pi*(x-a)/l);
    fpw = @(x) func_p_w(x).*cos((n-1)*pi*(x-a)/l);
    y(n)=integral(fw,a,b)/l;
    y_p(n)=integral(fpw,a,b)/l;
end

