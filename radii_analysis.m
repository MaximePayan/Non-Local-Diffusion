addpath functions
%% Data
close all

params = para;
params.delta = 4; % to test 1; 4; 2.445
if exist('intval','file')
    ipi = intval('pi');
else
    ipi = pi; 
end

params.domain = [0 2];
params.domain1 = [pi/4 pi/2];
params.domain2 = [pi/5 pi/2+0.25];

Ndiag = 5; % to test 10; 50; 100;
p = 1.7; % to test 0.5; 1; 2; 2.9
parameter_f = {@(x) x.^p,p,"pow"};

[B,C,q] = init_B(1,pi,Ndiag,params.domain,params.domain1,params.domain2);
parameter_B = {B,C,q};

%% without Intlab
disp("Precalculation without Intlab")
[B,C,q] = init_B(1,ipi,Ndiag-1,params.domain,params.domain1,params.domain2);
[RadiiBounds_prec,prec,M4,faraway_max_prec] = test_radii_anlysis_v2(params,pi,parameter_f,parameter_B);
[one_positivity,one_isolation]=spectrum_analysis(M4,RadiiBounds_prec,faraway_max_prec,true);
%% with Intlab
params_intval = params;
% Warning : If in params there are real numbers you want to be exact with
% (s.t. pi or sqrt(2)) you have to use them in intlab (ex:If params.a==sqrt(2),
% params_intval.a = intval('sqrt2') and not just intval(sqrt(2)).
if exist('intval','file')
    params_intval.domain = intval(params.domain);
    params_intval.a = intval(params.a);
    params_intval.b = intval(params.b);
    params_intval.c = intval(params.c);
    params_intval.d = intval(params.d);
    params_intval.delta = infsup(max(params.delta-dd,0),params.delta+dd);%intval(params.delta);
    params_intval.theta = intval(params.theta);
    params_intval.domain1 = [ipi/4 ipi/2];%intval(params.domain1);%
    params_intval.domain2 = [ipi/5 ipi/2+intval(0.25)];%intval(params.domain2);%
    [B_intval,C_intval,q_intval] = init_B(1,ipi,Ndiag,params_intval.domain,params_intval.domain1,params_intval.domain2);
    parameter_B_intval = {B_intval,C_intval,q_intval};
    parameter_f_intval = {@(x) x.^intval(p),intval(p),"pow"};
    disp("Rigorous calculations with Intlab")
    [RadiiBounds,prec,M4_intval,faraway_max] = test_radii_anlysis_v2test_radii_anlysis_v2(params_intval,ipi,parameter_f_intval,parameter_B_intval);
    [one_positivity,one_isolation,mu,d0]=spectrum_analysis(M4_intval,RadiiBounds,faraway_max,true);
end