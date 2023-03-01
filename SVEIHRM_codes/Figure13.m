clear all
load('plot_ori_delta')
vac=1;
sa=0;


Delta_table = readtable('Delta_1_1_100.csv');  % 
delta_variant=Delta_table{:,1};
vac_rate=1;
t = datetime(2020,2,15) + caldays(1:length(yyy));
tt=datenum(t);
hold on

% parameters
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
% 
% delta1=0.1;
% delta2=0.05;
% delta3=0.02;
% delta1_m=0.5;
% delta2_m=0.4;
% delta3_m=0.3;

parameters=[N;sigma;eta1;eta2;eta3;eta_m1;eta_m2;eta_m3;...
    gamma1;gamma2;gamma3;delta1;delta2;delta3;delta1_m;delta2_m;delta3_m];

ps=1;
k=1;
kk=1;
alpha1=400000;
alpha2=400000/4;
alpha3=400000/8;

for Mm=[0.9 1.0 1.1]
for i=1:100


para=[beta_ave*Mm tau_ave];
SEIRV_initial = yvalstrue(:,2);
tspan = [i+length(I)-1,i+length(I)];
if i==1
    soltrue = ode45(@(t,y)diffun_m(t,y,para,ps,parameters,alpha1,alpha2,alpha3),tspan,SEIRV_initial);
    yvalstrue_p = deval(soltrue,tspan);
    SEIRV_initial2=yvalstrue_p(:,2);
    yyy_m(1)=yvalstrue_p(4,2);
    yim_m(1)=yvalstrue_p(5,2);

else
    soltrue = ode45(@(t,y)diffun_m(t,y,para,ps,parameters,alpha1,alpha2,alpha3),tspan,SEIRV_initial2);
    yvalstrue_p = deval(soltrue,tspan);
    SEIRV_initial2=yvalstrue_p(:,2);
    yyy_m(i)=yvalstrue_p(4,2);
    yim_m(i)=yvalstrue_p(5,2);

end
sol_p(:,i)=yvalstrue_p(:,2);
PRtm(i)=beta_ave*tau_ave/mu*(sol_p(1,i)+delta1_m*sol_p(12,i)+delta2_m*sol_p(13,i))/N;
t=datetime(2020,2,8+tspan);
tt=datenum(t);
% subplot(1,2,2)
% ggt=plot(tt,yvalstrue_p(4,:)+yvalstrue_p(5,:),'c','LineWidth',3);


if Mm==0.9
    gg1_9=plot(tt,yvalstrue_p(5,:),'color',[1, 0.3 ,0.7],'LineWidth',2);
elseif Mm==1
    gg1_1=plot(tt,yvalstrue_p(5,:),'color',[1, 0 ,0],'LineWidth',2);
else
    gg1_11=plot(tt,yvalstrue_p(5,:),'color',[1, 0.7 ,0.3],'LineWidth',2);
end    
gg2=plot(tt,yvalstrue_p(4,:),'b','LineWidth',2);
plot([738521,738522],[4874,3917],'k','LineWidth',0.5)
gg3=plot(tt+1,delta_variant(i:i+1),'k','LineWidth',0.5);

    
hold on
end
real_covid(1)=4874;
real_covid(2:100)=delta_variant(1:99)';
app_covid=yyy_m+yim_m;
RL2error=sqrt((sum((real_covid-app_covid).^2))/sum(real_covid.^2))
% RMSEerror=sqrt(sum((real_covid-app_covid).^2)/length(real_covid))


end
xlim([738522 738622])
ylim([-10 8000])
legend([ gg1_9 gg1_1 gg1_11 gg2 gg3],'90% x Average R_t^m ','Average R_t^m','110% x Average R_t^m','The original virus','Actual')
title('Prediction of infectious according to transmission rate')
xlabel('100 days later') 
ylabel('The number of infectious')
datetick('x','mm/dd')
grid on


file_name='[Fig]Vaccine_predict_';
Mm_number=num2str(Mm);
extension='.png';
name=strcat(file_name,Mm_number,extension);
if sa==1
saveas(gcf,name)
end


