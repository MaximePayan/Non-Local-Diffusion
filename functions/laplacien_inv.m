function L = laplacien_inv(N,l,mat)
if exist('intval','file')
    ipi = intval('pi');
    Zeros = intval(zeros(N,1));
else
    ipi = pi;
    Zeros = zeros(N,1);
end

if nargin == 2 || mat == false
    L = Zeros;
    L(2:N) = -(ipi/l)^(-2)*(1:N-1)'.^(-2);
else
    L = Zeros*Zeros';
    L(2:N,2:N) = -(ipi/l)^(-2)*diag((1:N-1).^(-2));
end
end
