function [Y,Z1,Z2,r_min,r_max] = proof(a,b,sigma,eps,nu,u,gv,gpv,r,gpv_r,gppv_r,F,DF)
if exist('intval','file')
    ipi = intval('pi');
    a=intval(a);b=intval(b);
    eps=intval(eps); sigma=intval(sigma);
else
    ipi = pi;
end
l = b-a;
D = length(F)/2;
N_out = len(gv);
DK = 2*N_out;

%% Definition of A and W
% inverse of a truncation of DF

A_bar = inv([DF{1,1}(1:4*N_out-2,1:4*N_out-2),DF{1,2}(1:4*N_out-2,1:4*N_out-2)...
    ;DF{2,1}(1:4*N_out-2,1:4*N_out-2),DF{2,2}(1:4*N_out-2,1:4*N_out-2)]); %size 4(2N_out-1)x4(2N_out-1)
%truncation of this inverse
A_bar = A_bar([1:(2*N_out-1),2*(2*N_out-1)+(1:2*N_out-1)],...
    [1:2*N_out-1,2*(2*N_out-1)+(1:2*N_out-1)]); %size 2(2N_out-1)x2(2N_out-1)

Lap_inv = laplacien_inv(D,l,true);
[A,W] = approx_inverse(A_bar,u,gv,gpv,Lap_inv,eps,N_out,D);
%% Calcualation of the bound Y
Lap = laplacien(D,l,true);
A11Lap = A{1,1}*Lap;
normeA11Lap = norme_nu(A11Lap,nu,true);
A21Lap = A{2,1}*Lap;
normeA21Lap = norme_nu(A21Lap,nu,true);
normeu = norme_nu(u,nu);
Y = norme_nu([A{1,1}(:,1:2*N_out-1);A{1,2}(:,1:2*N_out-1)]*F(1:2*N_out-1),nu) + ...
    norme_nu([A{2,1}(:,1:2*N_out-1);A{2,2}(:,1:2*N_out-1)]*F(D+(1:2*N_out-1)),nu) ...
    + (normeA11Lap+normeA21Lap)*normeu*gv.Err;

if exist('intval','file')
    disp(['Y = ',mat2str(infsup(Y))])
else
    disp(['Y = ',num2str(Y)])
end

%% Calculation of the bound Z1

B11 = eye(D,DK) - (A{1,1}(:,1:DK+2*N_out-1)*DF{1,1}(1:DK+2*N_out-1,1:DK) ...
    + A{1,2}(:,1:DK+2*N_out-1)*DF{2,1}(1:DK+2*N_out-1,1:DK));
B12 = - (A{1,1}(:,1:DK+2*N_out-1)*DF{1,2}(1:DK+2*N_out-1,1:DK) ...
    + A{1,2}(:,1:DK+2*N_out-1)*DF{2,2}(1:DK+2*N_out-1,1:DK));
B21 = - (A{2,1}(:,1:DK+2*N_out-1)*DF{1,1}(1:DK+2*N_out-1,1:DK) ...
    + A{2,2}(:,1:DK+2*N_out-1)*DF{2,1}(1:DK+2*N_out-1,1:DK));
B22 = eye(D,DK) - (A{2,1}(:,1:DK+2*N_out-1)*DF{1,2}(1:DK+2*N_out-1,1:DK) ...
    + A{2,2}(:,1:DK+2*N_out-1)*DF{2,2}(1:DK+2*N_out-1,1:DK));

cZ1 = zeros(DK,1);
k = ksi(DK,nu);

if exist('intval','file')
    cZ1=intval(cZ1);
    k=intval(k);
end
% norm of B
for i=1:DK
    cZ1(i) = max(norme_nu(B11(:,i),nu)+norme_nu(B12(:,i),nu),...
        norme_nu(B21(:,i),nu) + norme_nu(B22(:,i),nu))/k(i);
end

Z1_finite = max(cZ1);

one = eye(2*N_out-1,1);

normeW = max(norme_nu(W{1,1},nu)+norme_nu(W{2,1},nu),...
    norme_nu(W{1,2},nu)+norme_nu(W{2,2},nu));
norme_errW = max(norme_nu(one-convo(W{1,1},gv.App,2*N_out-1),nu) + ...
    norme_nu(convo(W{2,1},gv.App,2*N_out-1),nu),...
    norme_nu(convo(W{1,1},convo(gpv.App,u,2*N_out-1),2*N_out-1)+eps.*W{1,2},nu)) + ...
    norme_nu(one-convo(W{2,1},convo(gpv.App,u,2*N_out-1),2*N_out-1)-eps.*W{2,2},nu);

Z1_tail = norme_errW + l^2/((DK-N_out+1)*ipi)^2*normeW*(sigma*norme_nu(one(1:N_out)-2.*u,nu)+1);

Z1 = Z1_finite + Z1_tail + (normeA11Lap+normeA21Lap)*(gv.Err+gpv.Err*normeu);

if exist('intval','file')
    disp(['Z1 = ',mat2str(infsup(Z1))])
else
    disp(['Z1 = ',num2str(Z1)])
end

%% Calculation on the bound Z2

Z2_a = (normeA11Lap+normeA21Lap)*(norm_err(gpv_r,nu) + gpv_r.Err);

Z2_b =(normeA11Lap + normeA21Lap)*(norme_nu(convo(gppv_r.App,u,2*N_out-1),nu)+...
    norm_err(gppv_r,nu)*r + (normeu+r)*gppv_r.Err);

normeA11 = norme_nu(A{1,1},nu,true);
normeA21 = norme_nu(A{2,1},nu,true);

Z2_c = 2*sigma*(norme_nu(A{1,1}(:,1:N_out)*u,nu)+norme_nu(A{2,1}(1:N_out)*u,nu)...
    + (normeA11 + normeA21)*r);

Z2 = max([Z2_a Z2_b Z2_c]);

if exist('intval','file')
    disp(['Z2 = ',mat2str(infsup(Z2))])
else
    disp(['Z2 = ',num2str(Z2)])
end
%% Check hypothesis 

if Z1 < 1
    disp 'Z1 < 1 : ok'
else
    disp 'Z1 >= 1 : not good'
end

test1 = 2*Y*Z2;
test2 = (1-Z1)^2;
if exist('intval','file')
    disp(['2Y*Z2 = ',mat2str(infsup(test1))])
    disp(['(1-Z1)^2 = ',mat2str(infsup(test2))])    
else
    disp(['2Y*Z2 = ',num2str(test1)])
    disp(['(1-Z1)^2 = ',num2str(test2)])
end

if test1 < test2
    disp '2Y*Z2 < (1-Z1)^2 : ok'
    r_min = (1-Z1 - sqrt(test2-test1))/Z2;
    if exist('intval','file')
        r_min = infsup(max(0,r_min.inf),r_min.sup);
    end
    r_max = (1-Z1)/Z2;
    if exist('intval','file')
        disp(['r_min = ',mat2str(infsup(r_min))])
        disp(['r_max = ',mat2str(infsup(r_max))])    
    else
        disp(['r_min = ',num2str(r_min)])
        disp(['r_max = ',num2str(r_max)])
    end
else
    disp '2Y*Z2 >= (1-Z1)^2 : not good'
    r_min = nan; r_max = nan;
end
end