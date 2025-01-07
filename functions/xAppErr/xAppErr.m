classdef xAppErr
    properties
        App 
        Err(1,1)
    end
    methods
        function x = xAppErr(app,err)
            if nargin == 1
                x.App = app;
            end
            if nargin == 2
                x.App = app;
                x.Err = err;
            end
        end
        function l = len(x)
            l = length(x.App);
        end
        function y = times(n,x)
            y = xAppErr;
            if isequal(size(n),[1,1])
                y.App = n.*x.App;
                y.Err = abs(n).*x.Err;
            else
                y.App = x.*n.App;
                y.Err = abs(x).*n.Err;
            end
        end
        function y = rdivide(x,n)
            y = xAppErr;
            y.App = x.App./n;
            y.Err = x.Err./abs(n);
        end
        function x = plus(x1,x2)
            x = xAppErr;
            l1 = len(x1); l2 = len(x2);
            if l1>l2
                x.App = x1.App;
                for i=1:l2
                    x.App(i) = x.App(i) + x2.App(i);
                end
            else
                x.App = x2.App;
                for i=1:l1
                    x.App(i) = x.App(i) + x1.App(i);
                end
            end
            x.Err = x1.Err + x2.Err;
        end
        function x = minus(x1,x2)
            x = xAppErr;
            l1 = len(x1); l2 = len(x2);
            if l1>l2
                x.App = x1.App;
                for i=1:l2
                    x.App(i) = x.App(i) - x2.App(i);
                end
            else
                x.App = -x2.App;
                for i=1:l1
                    x.App(i) = x.App(i) + x1.App(i);
                end
            end
            x.Err = x1.Err + x2.Err;
        end
    end
end