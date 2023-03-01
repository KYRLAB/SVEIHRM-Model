%%% Sensitivity analysis code for the vaccine.
for eff=[0.9 1. 1.1]
% eff=.9;
plt=1; %If 1, print the figure.
load('plot_ori_delta')
% 변수 정리
N=50000000;

sigma=1/4.1;
mu=1/4;

eta1=80/100;
eta2=18/100;
eta3=2/100;
eta_m1=80/100;
eta_m2=18/100;
eta_m3=2/100;

gamma1=1/14;
gamma2=1/23;
gamma3=1/30;


delta1=0.047;
delta2=0.002;
delta3=0.001;
delta1_m=0.5;
delta2_m=0.4;
delta3_m=0.3;


parameters=[N;sigma;eta1;eta2;eta3;eta_m1;eta_m2;eta_m3;...
    gamma1;gamma2;gamma3;delta1;delta2;delta3;delta1_m;delta2_m;delta3_m];
alpha1=alpha1*eff;
alpha2=alpha2*eff;
alpha3=alpha3*eff;

beta_betam(:,2)=beta_betam(:,2)*1.011;
sol=zeros(14,1);
ps=0;

for i=1:length(I)-8
    
    if vac==1
        ps=1;
    end
    para = optimvar('para',2,"LowerBound",0,"UpperBound",100);

    tspan=[i+7:i+8];
    SEIRV_initial=[50000000-40-9.93910598467827;40;0;9.93910598467827;0;0;0;0;0;0;0;0;0;0];
    True_interval_value=transpose(ytrue(i+7:i+8,1:2));
    if i==1
        myfcn = fcn2optimexpr(@ParaToODE,para,tspan,SEIRV_initial,ps,parameters,alpha1(i),alpha2(i),alpha3(i));
    else
        myfcn = fcn2optimexpr(@ParaToODE,para,tspan,SEIRV_initial2,ps,parameters,alpha1(i),alpha2(i),alpha3(i));
    end
    
    if i==1
        soltrue = ode45(@(t,y)diffun_m(t,y,beta_betam(i,:),...
            ps,parameters,alpha1(i),alpha2(i),alpha3(i)),tspan,SEIRV_initial);
        yvalstrue = deval(soltrue,tspan);
        SEIRV_initial2=yvalstrue(:,2);
        yyy(1)=yvalstrue(4,1);
    else
        soltrue = ode45(@(t,y)diffun_m(t,y,beta_betam(i,:),...
            ps,parameters,alpha1(i),alpha2(i),alpha3(i)),tspan,SEIRV_initial2);
        yvalstrue = deval(soltrue,tspan);
        SEIRV_initial2=yvalstrue(:,2);
        
        if i==300 %start the mutant
            SEIRV_initial2(3)=10;
            SEIRV_initial2(5)=41;
        end
    end
    
    yyy(i+1)=yvalstrue(4,2);
    yim(i+1)=yvalstrue(5,2);
    sol(:,i+1)=yvalstrue(:,2);
    t=datetime(2020,2,8+tspan);
    tt=datenum(t);
    if plt==1
%     plot(tt,II(tspan),'b')
    hold on
%     plot(tt,yvalstrue(4,:)+yvalstrue(5,:),'g','LineWidth',3)
%     plot(tt,yvalstrue(4,:),'r','LineWidth',3)
    hold on
    end
    
end


datetick('x','yy/mm/dd')


if plt==1
    grid on
    t = datetime(2020,2,15) + caldays(1:length(yyy));
    tt=datenum(t);
    plot(tt,yyy+yim,'LineWidth',2)
    xlim([datenum(datetime(2020,6,1)) 738521])

%     g2=plot(tt,yyy,'r','LineWidth',3)
%     g3=area(tt,ytrue(8:end,1:2),'FaceAlpha',.3,'EdgeAlpha',0)
%     colororder(['r';'g'])
%     g4=plot(tt,II(8:end),'b')
end
total=yyy+yim;
sum(total)
clear all
end
    
    legend({'90% x original vaccination dose','original vaccination dose','110% x original vaccination dose'},'location','northwest')
%     xlim([738300 738521])
%     legend([g1 g3 g4],{'Total number of infectious','Coronavirus','Matant virus','Actual'})
    title('Simulation of changes in past vaccination doses')
    xlabel('2020-06-01 ~ 2021-12-31') 
    ylabel('The number of infectious')
    datetick('x','mm/dd')
%     saveas(gcf,'[Fig]model_fitting2.png')
    
 