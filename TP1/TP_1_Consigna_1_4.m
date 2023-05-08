% Jorge Omar Sosa
% Sistema de Control II - 2023
% Trabajo Practico 1 - 1.4
% ===========================================================
% Definicón de las Matrices y los valores de Cada Variable
C= 680e-6; 
L= 23e-3; 
R= 11.6;  
A= [-R/L -1/L; 1/C 0];
B=[1/L; 0];
cT=[R 0];
% Tiempo de Integración:
t_I=1e-4;
% Tiempo de simulación
t_S=0.1;
% Pasos de la simulación
pasos=t_S/t_I;
t_inicial=0.05;
% Tomo solo valores mayores a 0,05 segundos
[~,lugar]=min(abs(t_inicial-tiempo));
il=[];
for ii=0:1:(length(tiempo)-1)
 if tiempo(ii+1)<=tiempo(lugar)
 il(ii+1)=0;
 else
 il(ii+1)=corriente_L(ii+1);
 end
end
% Matrices de los Datos
t=[pasos]; %tiempo
u=[pasos]; %entrada
corr_L=[pasos]; %corriente
%Condiciones Iniciales
X=[0;0]; 
y=[0]; 
% Datos Onda Cuadrada de Entrada
toc=50e-3; 
Ve=-1; 
for ii=0:1:pasos
 % Entrada del sistema
 u(ii+1)=Ve;
 if(tiempo(ii+1)<=toc)
 corr_L(ii+1)=0;
 else
 %variables de estados
 corr_L(ii+1)=X(1);
 %Sistema
 X_p=A*X+B*u(ii+1);
 X=X+t_I*X_p;
 Y=cT*X;
 end
 
end
figure(1);
% Sistema Original por Datos 
plot(tiempo,il,'red'); hold on
% Sistema Aproximado. 
plot(tiempo,corr_L,'blue');
title('Corriente L / Curva Aproximada')
xlabel('Tiempo [segundos]')
ylabel('Amplitud Corriente L')