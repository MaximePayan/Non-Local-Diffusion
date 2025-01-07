function [one_positivity,one_isolation,mu,d0]=spectrum_analysis(M,RadiiBounds,faraway_max,fig)

if exist('intval','file') && isintval(RadiiBounds(1))
    [Centers,Radii] = Gershgorin_disc(M);
    %Centers_intval is ordered, in decreasing order
    d_minus = real(Centers).inf - RadiiBounds.sup;
    d_plus = real(Centers).sup + RadiiBounds.sup;
    d_minus = [+Inf;d_minus];
    d_plus = [d_plus; faraway_max.sup];
else
    [Centers,Radii] = Gershgorin_disc(M);
    d_minus = real(Centers) - RadiiBounds;
    d_plus = real(Centers) + RadiiBounds;
    d_minus = [+Inf;d_minus];
    d_plus = [d_plus; faraway_max];
end
Isolated = [];
mu = nan;
for i=1:length(d_minus)-1
    if d_minus(i+1) > max(d_plus(i+1:end)) & d_plus(i) < min(d_minus(1:i))
        Isolated = [Isolated,i];
        if i == 1
            mu = max(d_plus(2:end));
            % if mu exists necessarily mu<d_minus(1) 
        end
    end
end

if nargout > 3
    d0 = Centers(1);
end
%% Tests
d0_minus = d_minus(2);
d0_plus = d_plus(1);
if isnan(mu)
    %disp('The first eigenvalue is not isolated')
    one_isolation = 0;
    j = 1;
    while not(ismember(1+j,Isolated)) && i+j < length(Centers)
        j = j+1;
    end
    if min(d_minus(1:1+j)) > 0 && max(d_plus(2+j:end)) < 0
        %disp("There are "+num2str(j)+" positive eigenvalues")
        one_positivity = 1;
    elseif max(d_plus(1:1+j)) < 0
        %disp("All eigenvalues are negative")
        one_positivity = 0;
    else
        %disp("We cannot conclude")
        one_positivity = nan;
    end
else
    one_isolation = 1;
    if mu<0 && d0_minus>0
        %disp("The first eigenvalue has a positive real part and is isolated")
        one_positivity = 1;
    elseif d0_plus < 0
        %disp("The first eigenvalue has a negative real part and is isolated")
        one_positivity = 0;
    else
        %disp("The first eigenvalue is isolated, but we cannot conclude on its sign")
        one_positivity = nan;
    end
end
%% Plots - centers with Radii Bounds
if nargin==4 && fig==true
theta = linspace(0,2*pi,1e2);
if exist('intval','file') && isintval(RadiiBounds(1))
    figure
    loglog(RadiiBounds.inf,'red','DisplayName','Bounds')
    hold on
    loglog(Radii.sup,'blue','DisplayName','Numerical Radii of diagonal part')
    title('Bounds with intlab')
    legend('Location','best')
    figure
    hold on
    for i=1:length(Centers)
        plot(real(Centers(i).mid),imag(Centers(i).mid),'b.');
        plot(Centers(i).mid + (RadiiBounds(i).sup+Centers(i).rad)*exp(1i*theta),'LineWidth',2,'Color','red')
    end
    plot([faraway_max.sup, faraway_max.sup],[-2,2],'LineWidth',2,'Color','black','DisplayName','Closest away circle')
    plot([0, 0],[-2,2],'b--','LineWidth',2,'DisplayName','x=0')
    title("Fisrt "+num2str(length(Centers))+" circles, with Intlab")
else
    figure
    loglog(RadiiBounds,'red','DisplayName','Bounds')
    hold on
    loglog(Radii,'blue','DisplayName','Numerical Radii')
    title('Bounds without intlab')
    legend('Location','best')

    tiledlayout(3,1)
    sizefont = 20;
    ax1=nexttile;
    hold on
    for i=1:length(Centers)
        plot(ax1,real(Centers(i)),imag(Centers(i)),'k.');
        plot(ax1,Centers(i) + RadiiBounds(i)*exp(1i*theta),'LineWidth',3,'Color',"#0072BD")
    end
    axis(ax1, "equal")
    limity = ylim(ax1);
    %ticks = xticks;
    plot(ax1,[0, 0],limity,'--','LineWidth',2,'DisplayName','x=0','Color',"#0072BD")
    %text(ax1,0-0.4*abs(ticks(2)-ticks(1)),0.9*limity(2),"y=0",'Color',"#0072BD",'FontSize',sizefont)
    plot(ax1,[faraway_max, faraway_max],limity,'--','LineWidth',2,'Color','black','DisplayName','x=M_{N,p}')
    plot(ax1,[mu, mu],limity,'LineWidth',2,'Color','black','DisplayName','x=\mu')
    %text(ax1,faraway_max_prec+0.1*abs(ticks(2)-ticks(1)),0.9*limity(2),"y=max_{i\geq 2N}(R_i+c_i)",'Color',"black",'FontSize',sizefont)
    ax1.FontSize = sizefont;

    ax2 = nexttile;
    hold on
    for i=1:5
        plot(ax2,real(Centers(i)),imag(Centers(i)),'k.');
        plot(ax2,Centers(i) + RadiiBounds(i)*exp(1i*theta),'LineWidth',3,'Color',"#0072BD")
    end
    plot(ax2,0,0,"-")
    xlim(ax2,[-3.5 1])
    axis(ax2, "equal")
    limity = ylim(ax2);
    %ticks = xticks;
    plot(ax2,[0, 0],limity,'--','LineWidth',2,'DisplayName','x=0','Color',"#0072BD")
    plot(ax2,[faraway_max, faraway_max],limity,'--','LineWidth',2,'Color','black','DisplayName','Closest away circle')
    plot(ax2,[mu, mu],limity,'LineWidth',2,'Color','black','DisplayName','x=\mu')
    %text(ax2,0-0.4*abs(ticks(2)-ticks(1)),0.9*limity(2),"y=0",'Color',"#0072BD",'FontSize',sizefont)
    ax2.FontSize = sizefont;

    ax3 = nexttile;
    hold on
    plot(ax3,real(Centers(1)),imag(Centers(1)),'k.');
    plot(ax3,Centers(1) + RadiiBounds(1)*exp(1i*theta),'LineWidth',3,'Color',"#0072BD")
    plot(ax3,0,0,"-")
    axis(ax3, "equal")
    limity = ylim(ax3);
    %ticks = xticks;
    plot(ax3,[0, 0],limity,'--','LineWidth',2,'DisplayName','x=0','Color',"#0072BD")
    %text(ax3,0+0.1*abs(ticks(2)-ticks(1)),0.9*limity(2),"y=0",'Color',"#0072BD",'FontSize',sizefont)
    ax3.FontSize = sizefont;
end
end
end
