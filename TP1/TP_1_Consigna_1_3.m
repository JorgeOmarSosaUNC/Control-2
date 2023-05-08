
% Jorge Omar Sosa
% Sistema de Control II - 2023
% Trabajo Practico 1 - 1.3
% ===========================================================
corriente_L = corriente_L(101:end);
tension_C = tension_C(101:end);
%Defino un nuevo tiempo para poder relacionarlo con las correspondientes
%tensiones y corrientes
aux=0:1e-4:((length(tension_C)-1)*1e-4);
tiempo_S=aux';
t_inicial=10e-3; % tiempo inicial en 10 mili Segundos
% busca el indice en el cual el tiempo inicial esta en la lista de valores del array "tiempo"
[~,punto]=min(abs(t_inicial-tiempo_S));
t_t1=tiempo_S(punto); % t1 el tiempo correspondiente al t_inicial
y_1=tension_C(punto); % y(t1)la tensión en C donde el tiempo=t_inicial
[~,punto]=min(abs(2*t_inicial-tiempo_S));
t_2t1=tiempo_S(punto); 
y_2=tension_C(punto); 
[~,punto]=min(abs(3*t_inicial-tiempo_S));
t_3t1=tiempo_S(punto); 
y_3=tension_C(punto);
 % Defino el Entrada Escalón
entr=stepDataOptions('InputOffset',0,'StepAmplitude',1);
% Considerando la salida como el último valor positivo de la Tensión en C:
K=tension_C(400)/entr.StepAmplitude;
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
plot(tiempo_S,tension_C,'red'); hold on
% Sistema Aproximado. 
step(G_s,entr) 
title('Curvas Medidas RLC / Curva Aproximada')
xlabel('Tiempo [segundos]')
ylabel('Amplitud Tensión en Capacitor')
damp(G_s)

