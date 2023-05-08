% Jorge Omar Sosa
% Sistema de Control II - 2023
% Trabajo Practico 1 - 2.3
% ===========================================================
tiempo = tiempo(104:15098);
VelocidadAng = VelocidadAng(104:15098); % valores positivos diferentes de 0
aux=tiempo;
tiempo_S=tiempo';
t_t1=tiempo_S(100); 
y_1=VelocidadAng(100); 
t_2t1=tiempo_S(5000); 
y_2=VelocidadAng(5000); 
t_3t1=tiempo_S(10000); 
y_3=VelocidadAng(10000); 
% Defino el Entrada Escalón del Sistema
entr=stepDataOptions('InputOffset',0,'StepAmplitude',1);
% último valor positivo de la Tensión en C:
K=VelocidadAng(end)/entr.StepAmplitude;
k1=(y_1/K)-1;
k2=(y_2/K)-1;
k3=(y_3/K)-1;
% Constantes para Chena
be=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3; 
alfa1=(k1*k2+k3-sqrt(be))/(2*(k1^2+k2)); 
alfa2=(k1*k2+k3+sqrt(be))/(2*(k1^2+k2)); 
beta=(k1+alfa2)/(alfa1-alfa2);
% Remplazo en los valores de las constantes de tiempo T1, T2 y T3:
T_1=-t_t1/log(alfa1); 
T_2=-t_t1/log(alfa2); 
T_3=beta*(T_1-T_2)+T_1; 
% Sistema Aproximado Final:
G_s=tf(K*[T_3 1],conv([T_2 1],[T_1 1]));
G_s=zpk(G_s);
% Sistema Original por Datos 
plot(tiempo_S,VelocidadAng,'r'); hold on
% Sistema Aproximado. 
step(G_s,entr) 
title('Curvas Medidas Motor 2023 / Curva Aproximada')
xlabel('Tiempo [segundos]')
ylabel('Velocidad Angular')

