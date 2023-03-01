clear all
load('plot_ori_delta')

figure
 t = datetime(2020,2,15) + caldays(1:length(yyy))+14;
    tt=datenum(t);
    plot(tt(1:length(Rt_real)),Rt_real,'b','LineWidth',1.5)
    hold on
    plot(tt(1:length(Rt_real)),Rtm,'r','LineWidth',1.5)
    hold on
    plot(tt(1:length(Rt_real)),ones(length(Rt_real)),'k','LineWidth',0.7)
    datetick('x','mm/dd')
    legend('R_t','R_t^m')
    xlim([738294 738529])
    ylim([-1. 3])

    title('R_t vs R_t^m')
    xlabel('2021-06-01 ~ 2021-12-31') 
    ylabel('Effective reproduction numbers')
    grid on