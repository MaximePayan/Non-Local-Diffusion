function [Fx,FDx] = handle_to_standard(f,fp,x)
%f and fp are two handle function, useful to return the handle function @(x)
%[f(x),fp(x)]
Fx = f(x);
FDx = fp(x);
end

