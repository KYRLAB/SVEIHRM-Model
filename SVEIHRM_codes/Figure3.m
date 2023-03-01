clear all


vact=100;
mutantt=200;
Pdate=400;

vac=0;
Rt=1.;
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


delta1=0.1;
delta2=0.05;
delta3=0.01;

delta1_m=0.5;
delta2_m=0.4;
delta3_m=0.3;


% delta1=0.047;
% delta2=0.002;
% delta1_m=0.5;
% delta2_m=0.4;

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


alpha1(vact+1:end)=N/300;
alpha2(vact+1:end)=N/1200;
alpha3(vact+1:end)=N/4800;



if vac==0
ps=0;
E1=40;Em1=0;I1=10;Im1=0;H11=0;H21=0;H31=0;R11=0;R21=0;R31=0;V11=0;V21=0;V31=0;
S1=N-E1-Em1-I1-Im1-H11-H21-H31-R11-R21-R31-V11-V21-V31;

initial_sol=[S1;E1;Em1;I1;Im1;H11;H21;H31;R11;R21;R31;V11;V21;V31];
sol=zeros(14,Pdate);

for Rt=[1.5 2 2.5 3 ]
k=1;
for i=1:Pdate
%%%%%%% No vaccine state. %%%%%%%%%%%%%%%
alpha1=0; 
alpha2=0;
alpha3=0;

delta1=0;
delta2=0;
delta3=0;
delta1_m=0;
delta2_m=0;
delta3_m=0;

parameters=[N;sigma;eta1;eta2;eta3;eta_m1;eta_m2;eta_m3;...
    gamma1;gamma2;gamma3;delta1;delta2;delta3;delta1_m;delta2_m;delta3_m];
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

% beta(i)=Rt*mu*N/(sol(1,i)+delta1*sol(12,i)+delta2*sol(13,i));
beta=para(1);
Rtt(i)=beta/mu/N*(sol(1,i)+delta1*sol(12,i)+delta2*sol(13,i)+delta3*sol(14,i));


end
t=1:Pdate+1;
plot(t,sol(4,:),'LineWidth',2)
% plot(Rtt,'LineWidth',2)
hold on
[Max_infectius, Max_day]=max(sol(4,:));
fprintf('Total number of infected : %d , Maximum infection rate : %f, Date : %d \n',sum(sol(4,:)),Max_infectius/7.245188e+06, Max_day)

end
legend('R_0=1.5','R_0=2.0','R_0=2.5','R_0=3.0')
title('Number of infectious individuals according to R_0')
xlabel('No vaccine or mutant virus for 400 days') 
ylabel('The number of infectious')
grid on

end






% 
% if vac==1
% ps=1;
% 
% E1=40;Em1=0;I1=10;Im1=0;H11=0;H21=0;H31=0;R11=0;R21=0;R31=0;V11=0;V21=0;
% S1=N-E1-Em1-I1-Im1-H11-H21-H31-R11-R21-R31-V11-V21;
% 
% initial_sol=[S1;E1;Em1;I1;Im1;H11;H21;H31;R11;R21;R31;V11;V21];
% sol=zeros(14,Pdate);
% 
% 
% k=1;
% for Rt=[3 ]
% for i=1:Pdate
% 
% sol(:,1)=initial_sol;
% 
% SEIRV_initial = sol(:,1);
% tspan = [i,i+1];
% 
% para=[Rt*mu Mm];%*N/(sol(1,530+i)+delta1*sol(12,503+i)+delta2*sol(13,530+i)) 1.]
% 
% if i==1
%     soltrue = ode45(@(t,y)diffun2(t,y,para,ps,delta1,delta2,delta1_m,delta2_m,...
%         alpha1(i),alpha2(i)),tspan,SEIRV_initial);
%     yvalstrue_p = deval(soltrue,tspan);
%     sol(:,i+1)=yvalstrue_p(:,2);
%     SEIRV_initial2=yvalstrue_p(:,2);
% else
%     soltrue = ode45(@(t,y)diffun2(t,y,para,ps,delta1,delta2,delta1_m,delta2_m,...
%         alpha1(i),alpha2(i)),tspan,SEIRV_initial2);
%     yvalstrue_p = deval(soltrue,tspan);
%     sol(:,i+1)=yvalstrue_p(:,2);
%     SEIRV_initial2=yvalstrue_p(:,2);
% end
% 
% 
% end
% t=1:Pdate+1;
% plot(t,sol(4,:),'LineWidth',2)
% hold on
% 
% end
% legend('R_t=1.5','R_t=2.0','R_t=2.5','R_t=3.0')
% title('Number of infectious individuals according to R_0')
% xlabel('400 days, No vaccine or mutant virus') 
% ylabel('The number of infectious')
% end
% grid on

% saveas(gcf,'paper_figure1.eps','epsc')



