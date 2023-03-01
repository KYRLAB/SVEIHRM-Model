function dydt = diffun_m(~,Equation,para,V_switch,parameters,alpha111,alpha222,alpha333)
dydt = zeros(14,1);

%%%% parameters
% parameters=[N;sigma;eta1;eta2;eta3;eta_m1;eta_m2;eta_m3;...
%     gamma1;gamma2;gamma3;delta1;delta2;delta3;delta1_m;delta2_m;delta3_m];
N=parameters(1);
sigma=parameters(2);
mu=1/4;

eta1=parameters(3);
eta2=parameters(4);
eta3=parameters(5);
eta_m1=parameters(6);
eta_m2=parameters(7);
eta_m3=parameters(8);

gamma1=parameters(9);
gamma2=parameters(10);
gamma3=parameters(11);


delta1=parameters(12);
delta2=parameters(13);
delta3=parameters(14);
delta1_m=parameters(15);
delta2_m=parameters(16);
delta3_m=parameters(17);

% Lambda=81;
% Mot=0.00006;
Lambda=788;
Mot=0.000114;

if V_switch==0
    alpha11=0;
    alpha22=0;
    alpha33=0;
end





S=Equation(1);
E=Equation(2);
Em=Equation(3);
I=Equation(4);
Im=Equation(5);
H1=Equation(6);
H2=Equation(7);
H3=Equation(8);
R1=Equation(9);
R2=Equation(10);
R3=Equation(11);
V1=Equation(12);
V2=Equation(13);
V3=Equation(14);

if S<400000
    S=0;
    alpha111=0;
end
if V1<400000
    V1=0;
    alpha222=0;
end
if V2<400000
    V2=0;
    alpha333=0;
end

SI_N=S*I/N;
SIm_N=S*Im/N;
V1I_N=V1*I/N;
V2I_N=V2*I/N;
V3I_N=V3*I/N;
V1Im_N=V1*Im/N;
V2Im_N=V2*Im/N;
V3Im_N=V3*Im/N;

if V_switch==0
    dydt(1) = Lambda-para(1)*SI_N-para(2)*para(1)*SIm_N-alpha11*S-Mot*S;
    dydt(2) = para(1)*SI_N-sigma*E+delta1*para(1)*V1I_N...
        +delta2*para(1)*V2I_N+delta3*para(1)*V3I_N-Mot*E;
    dydt(3) = para(2)*para(1)*SIm_N-sigma*Em+delta1_m*para(2)*para(1)*V1Im_N...
        +delta2_m*para(2)*para(1)*V2Im_N+delta3_m*para(2)*para(1)*V3Im_N-Mot*Em;
    dydt(4) = sigma*E-(eta1+eta2+eta3)*mu*I-Mot*I;
    dydt(5) = sigma*Em-(eta_m1+eta_m2+eta_m3)*mu*Im-Mot*Im;
    dydt(6) = eta1*mu*I+eta_m1*mu*Im-gamma1*H1-Mot*H1;
    dydt(7) = eta2*mu*I+eta_m2*mu*Im-gamma2*H2-Mot*H2;
    dydt(8) = eta3*mu*I+eta_m2*mu*Im-gamma3*H3-Mot*H3;
    dydt(9) = gamma1*H1-Mot*R1;
    dydt(10) = gamma2*H2-Mot*R2;
    dydt(11) = gamma3*H3-Mot*R3;
    dydt(12) = alpha11*S-alpha22*V1-delta1*para(1)*V1I_N-delta1_m*para(2)*para(1)*V1Im_N-Mot*V1;
    dydt(13) = alpha22*V1-alpha33*V2-delta2*para(1)*V2I_N-delta2_m*para(2)*para(1)*V2Im_N-Mot*V2;
    dydt(14) = alpha33*V2-delta3*para(1)*V3I_N-delta3_m*para(2)*para(1)*V3Im_N-Mot*V3;

else

    dydt(1) = -para(1)*SI_N-para(2)*para(1)*SIm_N-alpha111-Mot*S;
    dydt(2) = para(1)*SI_N-sigma*E+delta1*para(1)*V1I_N...
        +delta2*para(1)*V2I_N+delta3*para(1)*V3I_N-Mot*E;
    dydt(3) = para(2)*para(1)*SIm_N-sigma*Em+delta1_m*para(2)*para(1)*V1Im_N...
        +delta2_m*para(2)*para(1)*V2Im_N+delta3_m*para(2)*para(1)*V3Im_N-Mot*Em;
    dydt(4) = sigma*E-(eta1+eta2+eta3)*mu*I-Mot*I;
    dydt(5) = sigma*Em-(eta_m1+eta_m2+eta_m3)*mu*Im-Mot*Im;
    dydt(6) = eta1*mu*I+eta_m1*mu*Im-gamma1*H1-Mot*H1;
    dydt(7) = eta2*mu*I+eta_m2*mu*Im-gamma2*H2-Mot*H2;
    dydt(8) = eta3*mu*I+eta_m2*mu*Im-gamma3*H3-Mot*H3;
    dydt(9) = gamma1*H1-Mot*R1;
    dydt(10) = gamma2*H2-Mot*R2;
    dydt(11) = gamma3*H3-Mot*R3;
    dydt(12) = alpha111-delta1*para(1)*V1I_N-delta1_m*para(2)*para(1)*V1Im_N-alpha222-Mot*V1;
    dydt(13) = alpha222-delta2*para(1)*V2I_N-delta2_m*para(2)*para(1)*V2Im_N-alpha333-Mot*V2;
    dydt(14) = alpha333-delta3*para(1)*V3I_N-delta3_m*para(2)*para(1)*V3Im_N-Mot*V3;
end
