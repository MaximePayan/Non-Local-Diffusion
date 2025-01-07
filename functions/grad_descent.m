function [xopt,fopt,niter,gnorm,dx] = grad_descent(f,fp,x0,tol,maxiter,dxmin,alpha)
% grad_descent.m demonstrates how the gradient descent method can be used
% to solve a simple unconstrained optimization problem. Taking large step
% sizes can lead to algorithm instability. The variable alpha below
% specifies the fixed step size. Increasing alpha above 0.32 results in
% instability of the algorithm. An alternative approach would involve a
% variable step size determined through line search.
%
% This example was used originally for an optimization demonstration in ME
% 149, Engineering System Design Optimization, a graduate course taught at
% Tufts University in the Mechanical Engineering Department. A
% corresponding video is available at:
if nargin < 4
    % termination tolerance
    tol = 1e-6;
    % maximum number of allowed iterations
    maxiter = 1000;
    % minimum allowed perturbation
    dxmin = 1e-6;
    % step size ( 0.33 causes instability, 0.2 quite accurate)
    alpha = 0.1;
    % initialize gradient norm, optimization vector, iteration counter, perturbation
end
gnorm = inf; x = x0; niter = 0; dx = inf;
% define the objective function:
% plot objective function contours for visualization:
% redefine objective function syntax for use with optimization:
% gradient descent algorithm:
while and(gnorm>=tol, and(niter <= maxiter, dx >= dxmin))
    % calculate gradient:
    g = fp(x);
    gnorm = norm(g);
    % take step:
    xnew = x - alpha*g;
    % check step
    if ~isfinite(xnew)
        display(['Number of iterations: ' num2str(niter)])
        error('x is inf or NaN')
    end
    % update termination metrics
    niter = niter + 1;
    dx = norm(xnew-x);
    x = xnew; 
end
xopt = x;
fopt = f(xopt);
niter = niter - 1;
% define the gradient of the objective