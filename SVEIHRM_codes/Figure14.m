clear all

load('plot_omi')
% subplot(1,2,1)
grid on
    hold on
    t = datetime(2021,11,28) + caldays(1:length(yyy))-2;
    tt=datenum(t);
    g1=plot(tt,yyy+yim,'r','LineWidth',2); % total confirmed cases
    g2=plot(tt,yyy,'b','LineWidth',2);     % detla confirmed cases
    g3=area(tt,ytrue(1:end,1:2),'FaceAlpha',.3,'EdgeAlpha',0); % real
    colororder([0.,0.5,0.6;1.,0.5,0.5])
    g4=plot(tt,I_Cnt(1:end-1),'k');
    
    xlim([738487 738556])
    legend([g1 g3 g4],{'Total number of infectious','Delta variant','Omicron variant','Actual'},'location','northwest')
    title('Daily Delta and Omicron variant infections')
    xlabel('2021-11-28 ~ 2022-02-05') 
    ylabel('The number of infectious')
    datetick('x','mm/dd')
%     saveas(gcf,'[Fig]model_fitting2.png')
    
    
    figure
    gamma1=1/14;
    gamma2=1/23;
    gamma3=1/30;

    gamma=eta1*gamma1+eta2*gamma2+eta3*gamma3;
    Rt_C=beta*4;
    Rt_real=beta*4.*(sol(1,2:end)+delta1*sol(12,2:end)+delta2*sol(13,2:end))/N;
    Rtm=tau.*beta*4.*(sol(1,2:end)+delta1_m*sol(12,2:end)+delta2_m*sol(13,2:end))/N;
    %     +tau.*beta*4.*(sol(1,2:end)+delta1_m*sol(12,2:end)+delta2_m*sol(13,2:end))/N;
    t = datetime(2021,11,28) + caldays(1:length(yyy))-1;
    tt=datenum(t);
    plot(tt(1:length(Rt_real)),Rt_real,'b','LineWidth',2)
    hold on
    plot(tt(1:length(Rt_real)),Rtm,'r','LineWidth',2)
    hold on
    plot(tt(1:length(Rt_real)),(Rt_real.*(detect_delta(2:end-1)/100)'+Rtm.*(detect_omicron(2:end-1)/100)'),'LineWidth',2)
    hold on
    plot(tt(1:length(Rt_real)),ones(length(Rt_real)))
    datetick('x','mm/dd')
    legend('R_t^d','R_t^o')
    xlim([738516 738556])
    ylim([-0.5 3.5])

    title('R_t^d vs R_t^o')
    xlabel('2021-12-26 ~ 2022-02-05') 
    ylabel('Effective reproduction numbers')

    grid on
    
real_covid=II(1:end-1)+IIm(1:end-1);
app_covid=yyy+yim;
RL2error=sqrt(sum((real_covid-app_covid').^2)/sum(real_covid.^2))
RMSEerror=sqrt(sum((real_covid-app_covid').^2)/length(real_covid))