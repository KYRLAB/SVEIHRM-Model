%%% This is the data fitting code.


clear all
plt=1; %If 1, print the figure.
vt=434+28-7;
vac=1;
Mm=1;
% data import
CovidTable = readtable('covid_data_merge_20220108.csv','NumHeaderLines',3);  % skips the first three rows of data
VaccineData_t = readtable('vaccine.csv','NumHeaderLines',1);  % skips the first three rows of data
Vaccine1=VaccineData_t{:,5};
Vaccine2=VaccineData_t{:,6};
Vaccine3=VaccineData_t{:,7};


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

date_Cnt=CovidTable{:,1};
E_Cnt=CovidTable{:,2};
I_Cnt=CovidTable{:,6};
R_Cnt=CovidTable{:,4};
Im_Cnt=CovidTable{:,13};
PF=CovidTable{:,12};


E=zeros(length(E_Cnt)-1,1);
Em=zeros(length(E_Cnt)-1,1);
I=zeros(length(E_Cnt)-1,1);
Im=zeros(length(E_Cnt)-1,1);
H1=zeros(length(E_Cnt)-1,1);
H2=zeros(length(E_Cnt)-1,1);
H3=zeros(length(E_Cnt)-1,1);
R1=zeros(length(E_Cnt)-1,1);
R2=zeros(length(E_Cnt)-1,1);
R3=zeros(length(E_Cnt)-1,1);
V1=zeros(length(E_Cnt)-1,1);
V2=zeros(length(E_Cnt)-1,1);
V3=zeros(length(E_Cnt)-1,1);


alpha1=zeros(399+28+length(Vaccine1)-1,1);
alpha2=zeros(399+28+length(Vaccine2)-1,1);
alpha3=zeros(399+28+length(Vaccine3)-1,1);


alpha1(399+28:end)=Vaccine1;
alpha2(399+28:end)=Vaccine2;
alpha3(399+28:end)=0;


date=date_Cnt(2:end);

for i=1:length(E_Cnt)-1
    I(i,1)=I_Cnt(i+1)-I_Cnt(i);
    Im(i,1)=Im_Cnt(i+1);%-Im_Cnt(i);
    R1(i,1)=eta1*(R_Cnt(i+1)-R_Cnt(i));
    R2(i,1)=eta2*(R_Cnt(i+1)-R_Cnt(i));
    R3(i,1)=eta3*(R_Cnt(i+1)-R_Cnt(i));
end




pp=csaps(1:length(I),I,0.009);
pp2=csaps(1:length(I),Im,0.009);
ytrue=pp.coefs(:,4);
ytrue2=pp2.coefs(:,4);
% fnplt(pp)
% plot(ytrue)
% plot(I)
% hold on
ytrue(end+1)=sum(pp.coefs(end,:));
ytrue2(end+1)=sum(pp2.coefs(end,:));
ytrue(:,2)=ytrue2;
% ytrue(:,2)=ytrue.*PF(2:end)/100;
ytrue(:,3)=ytrue(:,1);
ytrue(:,1)=ytrue(:,1)-ytrue(:,2);
ytrue(ytrue<0)=0;


II=I;
% I=ytrue;
N=50000000;
S=N-E-Em-I-Im-H1-H2-H3-R1-R2-R3-V1-V2-V3;


sol=zeros(14,1);
ps=0;

for i=1:length(I)-8
    
    if vac==1
        ps=1;
    end
    para = optimvar('para',2,"LowerBound",0.01,"UpperBound",10);

    tspan=[i+7:i+8];
    SEIRV_initial=[50000000-40-9.93910598467827;40;0;9.93910598467827;0;0;0;0;0;0;0;0;0;0];
    True_interval_value=transpose(ytrue(i+7:i+8,1:2));
    if i==1
        myfcn = fcn2optimexpr(@ParaToODE,para,tspan,SEIRV_initial,ps,parameters,alpha1(i),alpha2(i),alpha3(i));
    else
        myfcn = fcn2optimexpr(@ParaToODE,para,tspan,SEIRV_initial2,ps,parameters,alpha1(i),alpha2(i),alpha3(i));
    end
    obj = sum(sum((myfcn - True_interval_value).^2));
    prob = optimproblem("Objective",obj);
    para0.para = [0.01 0];
    [parasol,sumsq] = solve(prob,para0);
    beta(i)=parasol.para(1);
    tau(i)=parasol.para(2);
    beta_betam(i,:)=[beta(i);tau(i)];
    error(i)=sumsq;
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
    g1=plot(tt,yyy+yim,'c','LineWidth',3);
    g2=plot(tt,yyy,'r','LineWidth',3);
    g3=area(tt,ytrue(8:end,1:2),'FaceAlpha',.3,'EdgeAlpha',0);
    colororder(['r';'g'])
    g4=plot(tt,II(8:end),'b');
    
    
    legend([g1 g3 g4],{'Total number of infectious','Coronavirus','Matant virus','Actual'})
    title('Daily coronavirus and mutant virus infections')
    xlabel('2020-02-16 ~ 2021-12-31') 
    ylabel('The number of infectious')
    datetick('x','mm/dd')
%     saveas(gcf,'[Fig]model_fitting2.png')
    
    
    figure
    gamma1=1/14;
    gamma2=1/23;
    gamma3=1/30;

    gamma=eta1*gamma1+eta2*gamma2+eta3*gamma3;
    Rt_real=beta*4.*(sol(1,2:end)+delta1*sol(12,2:end)+delta2*sol(13,2:end)+delta3*sol(14,2:end))/N;
    Rtm=tau.*beta*4.*(sol(1,2:end)+delta1_m*sol(12,2:end)+delta2_m*sol(13,2:end)+delta3_m*sol(14,2:end))/N;
    %     +tau.*beta*4.*(sol(1,2:end)+delta1_m*sol(12,2:end)+delta2_m*sol(13,2:end))/N;
    t = datetime(2020,2,15) + caldays(1:length(yyy))+14;
    tt=datenum(t);
    plot(tt(1:length(Rt_real)),Rt_real,'LineWidth',1.5)
    hold on
    plot(tt(1:length(Rt_real)),Rtm,'LineWidth',1.5)
    hold on
    plot(tt(1:length(Rt_real)),ones(length(Rt_real)))
    datetick('x','mm/dd')
    legend('R_t','R_t^m')
    title('R_t versus R_t^m')
    xlabel('2021-06-01 ~ 2021-12-31') 
    ylabel('Effective reproduction numbers')
    grid on
end


%%% Returns the mean Rt, beta, and tau over the last 4 weeks.

last=length(beta);
% figure
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
% save('beta_betum')
% save('plot_ori_delta')
