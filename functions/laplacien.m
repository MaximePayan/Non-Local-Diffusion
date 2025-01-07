function L=laplacien(N,l,mat)
if exist('intval','file') && isintval(l)
    ipi = intval('pi');
else
    ipi = pi;
end
if nargin == 2 || mat == false
    L = -(ipi/l)^2*(0:N-1)'.^2;
else
    L = -(ipi/l)^2*diag((0:N-1).^2);
end
end
