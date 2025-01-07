function y = Fourier_cos(w,x,l)
% w is a vector of fourier's coeff, x an abscisse, return the evaluation of the function
% represented by w in x
N = length(w);
y = w(1);
for k=1:N-1
    y = y + 2*w(k+1)*cos(k*pi*x/l);
end


