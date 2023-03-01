
clear all
for Rt=[1.4]
    for Mm=[1.3:0.1:1.8]
figure
vact=100;
mutantt=200; mswitch=1;
Pdate=700;
tvac=300;

vac=1;
mu=1/4;
beta=Rt*mu;

% Mm=1.7;

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





if vac==0 | vac==2
ps=0;
E1=40;Em1=0;I1=10;Im1=0;H11=0;H21=0;H31=0;R11=0;R21=0;R31=0;V11=0;V21=0;V31=0;
S1=N-E1-Em1-I1-Im1-H11-H21-H31-R11-R21-R31-V11-V21-V31;

initial_sol=[S1;E1;Em1;I1;Im1;H11;H21;H31;R11;R21;R31;V11;V21;V31];
sol=zeros(14,Pdate);


k=1;
for i=1:Pdate
alpha1=0;
alpha2=0;
alpha3=0;
delta1=0;
delta2=0;
delta3=0;
delta1_m=0;
delta2_m=0;
delta3_m=0;
sol(:,1)=initial_sol;

para=[beta Mm];

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
end


for vact=[ 100]
clearvars sol
alpha1=zeros(Pdate,1);
alpha2=zeros(Pdate,1);
alpha3=zeros(Pdate,1);

alpha1(vact+1:vact+300)=N/tvac;
alpha2(vact+1:end)=N/(tvac*4);
alpha3(vact+1:end)=N/(tvac*8);




if vac==1  | vac==2
ps=1;

E1=40;Em1=0;I1=10;Im1=0;H11=0;H21=0;H31=0;R11=0;R21=0;R31=0;V11=0;V21=0;V31=0;
S1=N-E1-Em1-I1-Im1-H11-H21-H31-R11-R21-R31-V11-V21-V31;

initial_sol=[S1;E1;Em1;I1;Im1;H11;H21;H31;R11;R21;R31;V11;V21;V31];
sol=zeros(14,Pdate);


k=1;

for i=1:Pdate

sol(:,1)=initial_sol;

SEIRV_initial = sol(:,1);
parameters=[N;sigma;eta1;eta2;eta3;eta_m1;eta_m2;eta_m3;...
    gamma1;gamma2;gamma3;delta1;delta2;delta3;delta1_m;delta2_m;delta3_m];
tspan = [i,i+1];

para=[beta Mm];

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
    if i==mutantt && mswitch==1
        
        sol(:,i+1)=yvalstrue_p(:,2);
        sol(3,i+1)=4;
        sol(5,i+1)=1;
        SEIRV_initial2=yvalstrue_p(:,2);
        SEIRV_initial2(3)=4;
        SEIRV_initial2(5)=1;

    end
end


end
t=1:Pdate+1;
subplot(1,2,1)
plot(t,sol(4,:)+sol(5,:),'LineWidth',2)
hold on

hold on
Rt_nor=beta*4.*(sol(1,1:end)+delta1*sol(12,1:end)+delta2*sol(13,1:end)+delta3*sol(14,1:end))/N;
Rt_mut=Mm.*beta*4.*(sol(1,1:end)+delta1_m*sol(12,1:end)+delta2_m*sol(13,1:end)+delta3*sol(14,1:end))/N;
end


grid on
fprintf('코로나 감염 비율 : %f (%%), 변이 바이러스 감염 비율 : %f (%%) \n',sum(sol(4,:))/N*100,sum(sol(5,:))/N*100)
end
title('Effect of Vaccine with mutant')
xlabel_name=sprintf('τ=%.1f',Mm);
xlabel(xlabel_name)
ylabel('The number of infectious people')
% limits=[0 6*10^5];
% ylim(limits)
% legend('vact=75','vact=100','vact=125')
% save_name=sprintf('vact=%d_Rt=%.1f_Mm=%.1f.png',vact,Rt,Mm)
% saveas(gcf,save_name)
% close
Rt_nor(Rt_nor<0)=0;
Rt_mut(Rt_mut<0)=0;
Rt_mut(1:mutantt-1)=0;
% figure
subplot(1,2,2)
plot(Rt_nor,linewidth=1.5)
hold on
plot(Rt_mut,linewidth=1.5)
limits=[0 2.2];
ylim(limits)
grid on
title('R_t verse R_t^m')
xlabel(xlabel_name) 
legend('R_t','R_t^m')
power_of_Nm=Rt_mut./Rt_nor;
Rt_nor(200)
Rt_mut(200)
power_of_Nm(200);
set(gcf,'position',[10,10,900,400])
% save_name=sprintf('Figure4_re_Mm=%.1f.eps',Mm)
% saveas(gcf,save_name,'epsc')
% % save_name=sprintf('vact=%d_Rt=%.1f_Mm=%.1f_%.2f_%.2f.eps',vact,Rt,Mm,sum(sol(4,:))/N*100,sum(sol(5,:))/N*100)
% % saveas(gcf,save_name,'epsc')
% close
    end
end
