clear all
load('plot_ori_delta')
fig1 = figure(1);




set(fig1, 'OuterPosition', [3, 270, 700, 420])

% subplot(1,2,1)
grid on
    t = datetime(2020,2,15) + caldays(1:length(yyy));
    tt=datenum(t);
    g1=plot(tt,yyy+yim,'r','LineWidth',2);
    hold on
    g2=plot(tt,yyy,'b','LineWidth',2);
    hold on

    g3=area(tt,ytrue(8:end,1:2),'FaceAlpha',.3,'EdgeAlpha',0);
    
    hold on
    colororder([0.,0.5,0.6;1.,0.5,0.5])
    g4=plot(tt,II(8:end),'k');
    grid on
    
    
    legend([g1 g3 g4],{'Total number of infected individuals','CoV','MuV (Delta variant)','Actual'},'Location','northwest')
    title('Daily coronavirus and mutant virus infections')
    xlabel('2020-02-16 ~ 2021-12-31') 
    ylabel('The number of infectious')
    datetick('x','mm/dd')
%%%%%%%%%%%% relative l2 error

real_covid=II(8:end);
app_covid=yyy+yim;
RL2error=sqrt(sum((real_covid-app_covid').^2)/sum(real_covid.^2))
RMSEerror=sqrt(sum((real_covid-app_covid').^2)/length(real_covid))
% RMSEPerror=RMSEerror/(sum(real_covid)/length(real_covid))*100


for i=0:4
    
    beta((last-6-(4-i)*7):last-(4-i)*7);
    Rtm((last-6-(4-i)*7):last-(4-i)*7);
    mean(Rtm((last-6-(4-i)*7):last-(4-i)*7))
    
    ave2(i+1)=sum(Rt_real((last-6-(4-i)*7):last-(4-i)*7))/7;
    beta_ave2(i+1)=sum(beta((last-6-(4-i)*7):last-(4-i)*7))/7;

end
% ave2=ave2*4;
ave2(1:4);
Rt_ave=mean(ave2(1:4));
tau_ave=mean(tau)
beta_ave=mean(beta_ave2(1:4));


