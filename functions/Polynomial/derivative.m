function Pp = derivative(P)
if not(class(P)=="Polynomial")
    P = Polynomial(P);
end
d = P.deg;
Pp = Polynomial;
if d>=1
    Pp.Value = (1:d).*P.Value(2:d+1);
end
end

