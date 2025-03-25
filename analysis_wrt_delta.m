addpath functions
Ddelta = linspace(0,4,1001);
% dd = 0.5*(4-0)/(length(Ddelta)-1);
V = []; params = para;
Ndiag = 50;

p=2;
parameter_f = {@(x) x.^p,p,"pow"};

if exist('intval','file')
    ipi = intval('pi');
else
    ipi = pi;
end
%% 
[B,C,q] = init_B(1,pi,Ndiag,params.domain,params.domain1,params.domain2);
parameter_B = {B,C,q};

if exist('intval','file')
    params_intval = para;
    params_intval.domain = intval(params.domain);
    params_intval.a = intval(params.a);
    params_intval.b = intval(params.b);
    params_intval.c = intval(params.c);
    params_intval.d = intval(params.d);
    params_intval.theta = intval(params.theta);
    params_intval.domain1 = [ipi/4 ipi/2];%intval(params.domain1);%
    params_intval.domain2 = [ipi/5 ipi/2+intval(0.25)];%intval(params.domain2);%
    [B_intval,C_intval,q_intval] = init_B(1,ipi,Ndiag,params_intval.domain,params_intval.domain1,params_intval.domain2);
    parameter_B_intval = {B_intval,C_intval,q_intval};
    parameter_f_intval = {@(x) x.^intval(p),intval(p),"pow"};
end

tic
for i=1:length(Ddelta)-1
    params.delta = Ddelta(i);
    delta = params.delta;
    [RadiiBounds,prec,M4,faraway_max] = test_radii_anlysis_v2(params,pi,parameter_f,parameter_B);
    if exist('intval','file')
        params_intval.delta = infsup(delta,Ddelta(i+1));
        delta = params_intval.delta;
        [RadiiBounds,prec,M4,faraway_max] = test_radii_anlysis_v2(params_intval,ipi,parameter_f_intval,parameter_B_intval,prec);
    end
    [one_positivity,one_isolation,mu,d0] = spectrum_analysis(M4,RadiiBounds,faraway_max);
    V = [V;delta, d0, RadiiBounds(1), mu, one_positivity, one_isolation];
end
toc
%%
color = ["#D95319","r","#77AC30","g","k"];
% dd = 0.5*(4-0)/(length(Ddelta)-1);
if exist('intval','file')
    save("first_eigv_wrt_delta_temp.mat",'Ddelta','V','params_intval','params','parameter_B','parameter_B_intval')
    color_index = V(:,5).mid + 2*V(:,6).mid + 1;
    color_index = fillmissing(color_index,'constant',5);
    %tiledlayout(1,2, TileSpacing="tight",Padding="loose")
    %ax1 = nexttile;
    ax1 = axes;
    sizefont = 20;
    %plot(real(V(:,1).mid),real(V(:,2).mid),'ko')
    hold on
    for i =1:length(Ddelta)-1
        rectangle('Position',[real(V(i,1).inf),real(V(i,2).inf-V(i,3).sup),V(i,1).sup-V(i,1).inf,real(V(i,2).sup-V(i,2).inf + 2*V(i,3).sup)],'EdgeColor',color(color_index(i)),'LineStyle','-','LineWidth',2)
    end
    xlim([0 4])
    %xlabel("\delta")
    %ax1.XTickLabel = {'','','',''};
    %ylabel("Re(d_0)","Rotation",0,"HorizontalAlignment","center")
    axis(ax1,'square')
    yax1 = ylim(ax1);
    ax1.FontSize = sizefont;
    ax1.Box = 'on';
    %ax1.XGrid = 'on';
    ax1.PositionConstraint = 'InnerPosition';
    ax1.InnerPosition = [0.15, 0.15, 0.75, 0.75];
    
    %ax2 = nexttile;
    ax2 = ax1;
    hold on
    for i =1:length(Ddelta)-1
        plot([V(i,1).inf,V(i,1).sup],[V(i,4).sup V(i,4).sup],'-','LineWidth',2,'Color',"#0072BD")
    end
    %xlim([0 4])
    %xlabel("\delta")
    %ylabel("\mu","Rotation",0,"HorizontalAlignment","left")
    %axis(ax2,'normal')
    %yax2 = ylim(ax2);
    %m_yax2 = (yax2(2)+yax2(1))/2;
    %demi_yax1 = (yax1(2)-yax1(1))/2;
    %ylim(ax2, [m_yax2-demi_yax1 m_yax2+demi_yax1]); 
    %ax2.FontSize = sizefont;
    %ax2.Box = 'on';
    %ax2.XGrid = 'on';
    %ax2.PositionConstraint = 'InnerPosition';
    %ax2.InnerPosition = [0.15, 0.15, 0.75, 0.75];

else
    save("first_eigv_wrt_delta_temp.mat",'Ddelta','V','params','parameter_B')
    color_index = V(:,4) + 2*V(:,5) + 1;
    color_index = fillmissing(color_index,'constant',5);
    tiledlayout(2,1)
    sizefont = 20;
    ax1 = nexttile;
    plot(ax1,real(V(:,1)),real(V(:,2)),'ko')
    xlabel("\delta")
    ylabel("Re(d_0)",Rotation=0)
    ax1.FontSize = sizefont;
    ax2 = nexttile;
    plot(ax2,real(V(:,1)),V(:,4),'bo')
    xlabel("\delta")
    ylabel("\mu",Rotation=0)
    ax2.FontSize = sizefont;
    warning("The graphs cannot be taken as the truth")
end