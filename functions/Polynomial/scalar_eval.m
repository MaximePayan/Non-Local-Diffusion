function y = scalar_eval(x,P)
if not(class(P)=="Polynomial")
    P = Polynomial(P);
end
if P.eq_to_zero
    y = 0;
else
%Polynomes are represented by the list of there coeff
%Pnum the numerator polynomial, Pden the denominator polynomial
    X = x'.^(0:deg(P));
    y = P.Value*X';
end
end