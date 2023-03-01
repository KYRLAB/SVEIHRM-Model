clear all

N=50000000;
x1=[50:10:120];
y1=[300 255 211 168 127 90 56 27];
x2=x1;
y2=[300 261 222 183 146 109 74 43];

       p1 =     0.01632 ; 
       p2 =      -6.678 ; 
       p3 =       595.8 ; 
       f1=@(x)  p1*x.^2 + p2*x + p3;
y11=f1(x1);
y111=N./y11;


       p1 =    0.008139 ; 
       p2 =      -5.058  ;
       p3 =       534.5 ;
       f2=@(x)  p1*x.^2 + p2*x + p3;
y22=f2(x1);
y222=N./y22;

fig1 = figure(1);
set(fig1, 'OuterPosition', [3, 270, 840, 420])

subplot(1,2,1)
scatter(x1,y1)
hold on
plot(x1,y11)
xlabel('VS')
ylabel('1/VD')
title('1/VD according to VS with R_0=1.5')
text(75,193,'\leftarrow 0.01632 x^2 - 6.678 x + 595.8')
legend('','Curve Fitting with Quadratic Function')
grid on

subplot(1,2,2)
scatter(x2,y2)
hold on 
plot(x2,y22)
xlabel('VS')
ylabel('1/VD')
title('1/VD according to VS with R_0=1.8')
text(75,205,'\leftarrow 0.008139 x^2 - 5.058 x + 534.5')
legend('','Curve Fitting with Quadratic Function')
grid on

% save_name='VD according to VS'
% saveas(gcf,save_name,'epsc')
