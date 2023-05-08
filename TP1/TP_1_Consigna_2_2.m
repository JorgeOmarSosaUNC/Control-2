% Jorge Omar Sosa
% Sistema de Control II - 2023
% Trabajo Practico 1 - 2.2
% ===========================================================
t_S=1e-7; 
tF=0.04; 
u=12; 
TL=2.128e-5; %Torque
jj=0; 
X=-[0; 0; 0]; % Vector de Omega, Wr y ia
ii=0; 
for t=0:t_S:tF
ii=ii+1;
X=modmotor2(t_S, X, u, TL);
x1(ii)=X(1); % Omega
x3(ii)=X(3); % ia
ent(ii)=u;
end
t=0:t_S:tF;
subplot(3,1,1);hold on;
plot(t,x1,'r');title('Salida y, \omega_t');
subplot(3,1,2);hold on;
plot(t,x3,'b');title('Corriente I_a');
subplot(3,1,3);hold on;
plot(t,ent,'g');title('Entrada V_a');
xlabel('Tiempo [Seg.]');
