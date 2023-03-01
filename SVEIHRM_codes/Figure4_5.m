clear all
% figure
% for Rt=1.5:0.1:2.5
for Rt=2
vact=100;
mutantt=200;
Pdate=250;
tvac=300;

vac=2;
% Rt=1.6;%% 1.2744   1.4018   1.1470
mu=1/4;
Mm=1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters
N=50000000;

sigma=1/4.1;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E=zeros(Pdate,1);
Em=zeros(Pdate,1);
I=zeros(Pdate,1);
Im=zeros(Pdate,1);
H1=zeros(Pdate,1);
H2=zeros(Pdate,1);
H3=zeros(Pdate,1);
R1=zeros(Pdate,1);
R2=zeros(Pdate,1);
R3=zeros(Pdate,1);
V1=zeros(Pdate,1);
V2=zeros(Pdate,1);
V3=zeros(Pdate,1);


alpha1=zeros(Pdate,1);
alpha2=zeros(Pdate,1);
alpha3=zeros(Pdate,1);


alpha1(vact+1:end)=N/tvac;
alpha2(vact+1:end)=N/(tvac*4);
alpha3(vact+1:end)=N/(tvac*8);


if vac==0 | vac==2
ps=0;
E1=40;Em1=0;I1=10;Im1=0;H11=0;H21=0;H31=0;R11=0;R21=0;R31=0;V11=0;V21=0;V31=0;
S1=N-E1-Em1-I1-Im1-H11-H21-H31-R11-R21-R31-V11-V21-V31;

initial_sol=[S1;E1;Em1;I1;Im1;H11;H21;H31;R11;R21;R31;V11;V21;V31];
sol=zeros(14,Pdate);


k=1;
for i=1:Pdate

sol(:,1)=initial_sol;

para=[Rt*mu Mm];

SEIRV_initial = sol(:,1);
tspan = [i,i+1];

if i==1
    soltrue = ode45(@(t,y)diffun_m(t,y,para,ps,parameters,alpha1,alpha2,alpha3),tspan,SEIRV_initial);
    yvalstrue_p = deval(soltrue,tspan);
    sol(:,i+1)=yvalstrue_p(:,2);
    SEIRV_initial2=yvalstrue_p(:,2);
else
    soltrue = ode45(@(t,y)diffun_m(t,y,para,ps,parameters,alpha1,alpha2,alpha3),tspan,SEIRV_initial2);
    yvalstrue_p = deval(soltrue,tspan);
    sol(:,i+1)=yvalstrue_p(:,2);
    SEIRV_initial2=yvalstrue_p(:,2);
end

t=1:Pdate+1;


end

plot(t,sol(4,:),'LineWidth',2)

hold on
criterion=sum(sol(4,:));
end
% legend('R_t=1.5','R_t=2.0','R_t=2.5','R_t=3.0')
% title('No Vaccine and No mutant')
% xlabel('400 days') 
% ylabel('The number of infectious')


clearvars sol

alpha1=zeros(Pdate,1);
alpha2=zeros(Pdate,1);
alpha3=zeros(Pdate,1);




if vac==1  | vac==2
ps=1;

E1=40;Em1=0;I1=10;Im1=0;H11=0;H21=0;H31=0;R11=0;R21=0;R31=0;V11=0;V21=0;V31=0;
S1=N-E1-Em1-I1-Im1-H11-H21-H31-R11-R21-R31-V11-V21-V31;

initial_sol=[S1;E1;Em1;I1;Im1;H11;H21;H31;R11;R21;R31;V11;V21;V31];
sol=zeros(14,Pdate);


k=1;
for tvac=[300 200 100]
    alpha1(vact+1:end)=N/tvac;
    alpha2(vact+1:end)=N/(tvac*4);
    alpha3(vact+1:end)=N/(tvac*8);

parameters=[N;sigma;eta1;eta2;eta3;eta_m1;eta_m2;eta_m3;...
    gamma1;gamma2;gamma3;delta1;delta2;delta3;delta1_m;delta2_m;delta3_m];
for i=1:Pdate

sol(:,1)=initial_sol;

SEIRV_initial = sol(:,1);
tspan = [i,i+1];

para=[Rt*mu Mm];


if i==1
    soltrue = ode45(@(t,y)diffun_m(t,y,para,ps,parameters,alpha1(i),alpha2(i),alpha3(i)),tspan,SEIRV_initial);
    yvalstrue_p = deval(soltrue,tspan);
    sol(:,i+1)=yvalstrue_p(:,2);
    SEIRV_initial2=yvalstrue_p(:,2);
else
    soltrue = ode45(@(t,y)diffun_m(t,y,para,ps,parameters,alpha1(i),alpha2(i),alpha3(i)),tspan,SEIRV_initial2);
    yvalstrue_p = deval(soltrue,tspan);
    sol(:,i+1)=yvalstrue_p(:,2);
    SEIRV_initial2=yvalstrue_p(:,2);
end


end
t=1:Pdate+1;
plot(t,sol(4,:),'LineWidth',2)
hold on
fprintf('Ratio to non-vaccinated status : %f (%%) \n',sum(sol(4,:))/criterion*100)

end
% legend('R_t=1.5','R_t=2.0','R_t=2.5','R_t=3.0')
% title('Vaccine and No mutant')
% xlabel('400 days') 
% ylabel('The number of infectious')

end

text(142,3392420,'\leftarrow Criterion(100%)')
text(140,2773070,'\leftarrow 77.37%')
text(139,2493160,'\leftarrow 67.94%')
text(140,1602080,'\leftarrow 46.66%')
legend('No vaccine','VD=1/300','VD=1/200','VD=1/100')
title('Effectiveness of vaccines according to VD without mutant viruses')
xlabel('For 250 days when R_0=2.0') 
ylabel('The number of infectious people')
% xlabel(['R_0=',num2str(Rt)])
% limits=[0 6*10^6];
% ylim(limits)
grid on
% save_name=sprintf('VD=%d_Rt=%.1f_Mm=%.1f.eps',vact,Rt,Mm)
% saveas(gcf,save_name,'epsc')
% close
end

% 7/21 260명 (특이사항)