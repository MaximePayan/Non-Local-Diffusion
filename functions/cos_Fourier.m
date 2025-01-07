function w = cos_Fourier(f,N,a,b)
%COS_FOURIER Summary of this function goes here
%   Detailed explanation goes here
w = zeros(N,1);
l=b-a;
w(1) = integral(f,a,b)/l;
for n=2:N
    fn = @(x) f(x).*cos((n-1)*pi*(x-a)/l);
    w(n)=integral(fn,a,b)/l;
end


end

