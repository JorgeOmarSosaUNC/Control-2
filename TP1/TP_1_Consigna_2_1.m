% Jorge Omar Sosa
% Sistema de Control II - 2023
% Trabajo Practico 1 - 2.1
% ===========================================================
Laa=366e-6;
J=5e-9;
Ra=55.6;
B=0;
Ki=6.49e-3;
Km=6.53e-3;
den=[Laa*J Ra*J+Laa*B Ra*B+Ki*Km ];
num=Ki;
wr_va=tf(num,den);
sys=zpk(wr_va);
% Simulación
t_S=1e-7; %tiempo de simulación
tF=0.04; 
tF1=3e-7; 
u=12; 
TL=1e-5; %Torque inicial 
TLfin=1e-4; %Torque final 
TLi=1e-15; %Paso de cada Torque
jj=0;
while TL<TLfin
ii=0; 
X=-[0; 0]; % Vector de Omega y Wr
TL=TL+jj*TLi;
for t=0:t_S:tF1
ii=ii+1;
X=modmotor1(t_S, X, u, TL);
x1(ii)=X(1); %Omega
end
if X(1)<0
TL=TL-jj*TLi;
break;
end
jj=jj+1;
end
% Simulacion para grafico
ii=0;X=-[0; 0];
for t=0:t_S:tF
ii=ii+1;
X=modmotor1(t_S, X, u, TL);
x1(ii)=X(1);
end
tvector = linspace(0,tF,ii);
figure(1)
plot(tvector,x1);