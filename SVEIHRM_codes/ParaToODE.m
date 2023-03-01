function solparameter = ParaToODE(para,tspan,SEIRV_initial,V_switch,parameters,alpha1,alpha2,alpha3)
sol = ode45(@(t,Equation)diffun_m(t,Equation,para,V_switch,parameters,alpha1,alpha2,alpha3),tspan,SEIRV_initial);
solparameter = deval(sol,tspan);
solparameter=solparameter(4:5,:);
end