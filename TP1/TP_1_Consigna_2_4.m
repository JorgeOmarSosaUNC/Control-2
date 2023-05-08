% Jorge Omar Sosa
% Sistema de Control II - 2023
% Trabajo Practico 1 - 2.4
% ===========================================================
%Constantes del PID 
Kp=2;
Ki=100;
Kd=0;
X=-[0; 0; 0; 0];ii=0;t_etapa=1e-7;titaRef=1;tF=.2;
Ts=t_etapa; 
A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts); 
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts); 
C1=Kd/Ts; 
e=zeros(round(tF/t_etapa),1); 
u=0;u_max=12; 
for t=0:t_etapa:tF 
 ii=ii+1;k=ii+2; 
 X=modmotor4(t_etapa, X, u); 
 e(k)=titaRef-X(4); %ERROR 
 u=u+A1*e(k)+B1*e(k-1)+C1*e(k-2); %PID 
 x1(ii)=X(1);% Omega 
 x2(ii)=X(2);% wp 
 x3(ii)=X(3);% ia 
 x4(ii)=X(4);% tita 
 u=u_max*tanh(u/u_max); % termina la acción de control. 
 acc(ii)=u; 
end
 
t=0:t_etapa:tF; 
subplot(3,1,1);hold on; 
plot(t,x4);title('Salida posición \tita_t'); 
subplot(3,1,2);hold on; 
plot(t,x3);title('Corriente I_a'); 
subplot(3,1,3);hold on; 
plot(t,acc);title('Entrada V_a'); 
xlabel('Tiempo [Seg.]');