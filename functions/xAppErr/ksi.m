function k = ksi(N,nu)
k = ones(N,1);
for i=2:N
    k(i)=2*nu^(i-1);
end
end