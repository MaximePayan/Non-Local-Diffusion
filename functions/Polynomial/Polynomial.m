classdef Polynomial
    properties
        Value double {mustBeNumeric} = 0
    end
    methods
        function P = Polynomial(val)
            if nargin ==1 
                P.Value = val;
            end
        end
        function veracity = eq_to_zero(P)
            veracity = true;
            for i=1:length(P.Value)
                if P.Value(i) ~= 0
                    veracity = false;
                    break
                end
            end
        end
        function d = deg(P)
            if P.eq_to_zero
                d = -Inf;
            else
                d=0;
                while P.Value(end-d) == 0 && d>=0
                    d = d + 1;
                end
                d = length(P.Value)-d-1 ;
            end
        end
        function R = plus(P,Q)
            d1 = deg(P); d2=deg(Q);
            if d1 > d2
                R = P;
                for i=0:d2
                    R.Value(i+1) = R.Value(i+1) + Q.Value(i+1);
                end
            else
                R = Q;
                for i=0:d1 
                    R.Value(i+1) = R.Value(i+1) + P.Value(i+1);
                end
            end
        end
        function R = times(a,P)
            R = Polynomial;
            if isfloat(a)
                R.Value= a.*P.Value;
            else
                R.Value = P.*a.Value;
            end

        end
        function R = minus(P,Q)
            d1 = deg(P); d2=deg(Q);
            if d1 > d2
                R = P;
                for i=0:d2
                    R.Value(i+1) = P.Value(i+1) - Q.Value(i+1);
                end
            else
                R = (-1).*Q;
                for i=0:d1 
                    R.Value(i+1) = P.Value(i+1) - Q.Value(i+1);
                end
            end
        end
        function R = mtimes(P,Q)
            R = Polynomial;
            R.Value = conv(P.Value,Q.Value);
        end
        function veracity = eq(P,Q)
            veracity = eq_to_zero(P-Q);
        end
        function R = abs(P)
            R = Polynomial(abs(P.Value));
        end
    end
end