addpath functions
clear all
load first_eigv_wrt_delta_temp.mat

Ndiag = size(parameter_B{1},1);
D = 50;
%dd = 0.5*(Ddelta(1)-Ddelta(end))/(length(Ddelta)-1);

lambda = [];
lambda_p = [];
bound_lambda_p = [];

alpha = 0.01; %norm's parameter

if exist('intval','file')
    ipi = intval('pi');
    B = parameter_B_intval{1};
    C = parameter_B_intval{2};
    q = parameter_B_intval{3};
    index_nan = find(ismissing((V(:,5).mid)));
else
    ipi = pi;
    B = parameter_B{1};
    C = parameter_B{2};
    q = parameter_B{3};
    parameter_B_intval = parameter_B;
    index_nan = find(abs(V(:,5)) >= 0);
end

for i=1:length(index_nan)
    delta = V(index_nan(i),1);
    if exist('intval','file')
        params.delta = delta.mid;
        params_intval.delta = delta;
    else
        params.delta = delta;
        params_intval = params;
    end
    [~,~,~,~,prec] = script_nonlocal_diff_v2(parameter_B,params,D,alpha);
    [A,E_maj,chi_q,bounds,~] = script_nonlocal_diff_v2(parameter_B_intval,params_intval,D,alpha,prec);
    A_app = prec{3}; Utilde = prec{2}; 
    Y = bounds(1).sup
    Z1 = bounds(2).sup
    Z2 = bounds(3).sup
    r_min = bounds(4).sup
    r_max = bounds(5).inf
    if exist('intval','file') && r_min >=0
        u = infsup(Utilde(2:Ndiag+1)-r_min,Utilde(2:Ndiag+1)+r_min);
        d0 = midrad(Utilde(1),r_min);
    else
        u = Utilde(2:Ndiag+1);
        d0 = Utilde(1);
    end
    lp = -A(1,Ndiag+D+2:2*Ndiag+D+1)*B*u;
    if Z1 + Z2*r_min < 1
        b1 = C*abs(A(1,Ndiag+D+2:2*Ndiag+D+1))*chi_q*r_min + norme_alpha(abs(A(2:Ndiag+1,Ndiag+D+2:2*Ndiag+D+1))*chi_q,alpha)*r_min;
        b2 = norme_alpha(abs(A(Ndiag+D+2:2*Ndiag+D+1,Ndiag+D+2:2*Ndiag+D+1))*chi_q,alpha)*r_min+2*C*E_maj*r_min;
        b3 = 2*C*E_maj*abs(u)'*chi_q*r_min + (Z1+Z2*r_min)*norme_RL1L1(A_app*[0;zeros(Ndiag,1);B*u],alpha); 
        bound_lp = (b1+b2+b3)/(1-(Z1+Z2*r_min));
    else
        bound_lp = nan;
    end
    lambda = [lambda, d0];
    lambda_p = [lambda_p, lp];
    bound_lambda_p = [bound_lambda_p,bound_lp];
end
%%
fig  = true;
if fig
    if exist('intval','file')
        figure
        sizefont = 20;
        ax1 = axes;
        hold on
        for i =1:length(index_nan)
            delta = V(index_nan(i),1);
            rectangle(ax1,'Position',[delta.inf, ...
                lambda_p(i).inf-bound_lambda_p(i).sup,2.*delta.rad,2.*(lambda_p(i).rad+bound_lambda_p(i).sup)],'EdgeColor','b','LineStyle','-','LineWidth',2)
        end
        
        xlabel('\delta',FontSize=33)
        ylabel("Re(d_0')",'Rotation',0,FontSize=33,HorizontalAlignment='right')
        ax1.FontSize = sizefont;
        xlim(ax1,[V(index_nan(1),1).inf,V(index_nan(end),1).sup])
        axis(ax1,'normal')
        ax1.Box = 'on';
        ax1.PositionConstraint = 'InnerPosition';
        ax1.InnerPosition = [0.15, 0.15, 0.75, 0.75];
        figure
        ax2 = axes;
        hold on
        for i=1:length(index_nan)
            delta = V(index_nan(i),1);
%             p1 = lambda_p(i).inf-bound_lambda_p(i).sup;
%             p2 = lambda_p(i).sup+bound_lambda_p(i).sup;
%             plot(ax2,delta.mid,lambda(i).mid, 'bo' )
%             draw_circ_arc([delta.mid,lambda(i).mid],0.001,[delta.mid,lambda(i).mid]+delta.rad*[1,p1]./sqrt(1+p1^2),...
%                 [delta.mid,lambda(i).mid]+delta.rad.*[1,p2]./sqrt(1+p2^2));
            rectangle(ax2,'Position',[delta.inf, ...
                lambda(i).inf,2.*delta.rad,2.*lambda(i).rad],'EdgeColor','b','LineStyle','-','LineWidth',2)
        end
        xlabel('\delta',FontSize=33)
        ylabel("Re(d_0)",'Rotation',0,FontSize=33,HorizontalAlignment='right')
        ax2.FontSize = sizefont;
        xlim(ax2,[V(index_nan(1),1).inf,V(index_nan(end),1).sup])
        axis(ax2,'normal')
        ax2.Box = 'on';
        ax2.PositionConstraint = 'InnerPosition';
        ax2.InnerPosition = [0.15, 0.15, 0.75, 0.75];

    else
        tiledlayout(2,1)
        sizefont = 20;
        ax1 = nexttile;
        plot(ax1,V(index_nan,1),lambda_p-bound_lambda_p,'b-',DisplayName='inf')
        hold on
        plot(ax1,V(index_nan,1),lambda_p,'b--',DisplayName='mid')
        plot(ax1,V(index_nan,1),lambda_p+bound_lambda_p,'b-',DisplayName='sup')
        xlabel('\delta',FontSize=33)
        ylabel("Re(d_0')",'Rotation',0,FontSize=33,HorizontalAlignment='right')
        ax1.FontSize = sizefont;
        
        ax2 = nexttile;
        plot(ax2,V(index_nan,1),lambda,'ko')
        xlabel('\delta',FontSize=33)
        ylabel("Re(d_0)",'Rotation',0,FontSize=33,HorizontalAlignment='right')
        ax2.FontSize = sizefont;
        warning("The graphs cannot be taken as the truth")
    end
end