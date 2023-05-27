%Codigo del motor de CC
% Trabajo Practico 2 - Item 1 - 2023
% Alumno: Jorge Omar Sosa
clear all, close all, clc

%Variables
Laa= 5e-3;
J= 0.004;
Ra= 0.2;
B= 0.005;
Ki= 6.5e-5;
Km= 0.055;

%Torque (Dato)
T1=1.15e-3; %Va a ir variando para pi/2 y -pi/2


%Matrices (Luego de resolver usando variables de estado) 
Mat_A=[-Ra/Laa -Km/Laa 0; Ki/J -B/J 0; 0 1 0];
Mat_B=[1/Laa; 0; 0];
Mat_C=[0 0 1];% x1=ia   x2=wr   x3=theta
Mat_D=[0];

%Matrices ampliadas
Mat_Aa=[Mat_A zeros(3,1); -Mat_C 0];
Mat_Ba=[Mat_B; 0];
Mat_Cc=[Mat_C 0];

%Diseño del controlador por LQR ver si esto cambia el tiempo de cambio
Q=diag([1 1/10000 100 10000000000000]);
R=10.0; %Si achico R tiende a ser mas rapida la respuesta, o tambien voy aumentando Q

%Construccion del Hamiltoniano para el calculo del controlador
Ha=[Mat_Aa -Mat_Ba*inv(R)*Mat_Ba'; -Q -Mat_Aa'];
[n,va]=size(Ha);

[V,D]=eig(Ha);MX1X2=[];
for ii=1:n
    if real(D(ii,ii))<0
        MX1X2=[MX1X2 V(:,ii)];
    end
end

MX1=MX1X2(1:n/2,:);
MX2=MX1X2(n/2+1:end,:);

P=real(MX2*inv(MX1));%Tomo la parte real por un tema de residuos
Ka=inv(R)*Mat_Ba'*P;
K=Ka(1:3); 
KI=-Ka(4);
eig(Mat_Aa-Mat_Ba*Ka)%Verifico polos con parte real negativa
%Fin cálculo del controlador

%Para la simulación
delta_t=1e-4; 

tiempo=10;
pasos=round(tiempo/delta_t);
Ci=[0 0 0 0];
t=0:delta_t:(tiempo-delta_t);

%Referencia
ref=(pi/2)*square(2*pi*t/4);

T11=(T1/2)+(T1/2)*square(2*pi*t/4 );%Función para ir variando el torque

%Condiciones iniciales
x=zeros(4,pasos);
x(1,1)=Ci(1);%ia
x(2,1)=Ci(2);%w
x(3,1)=Ci(3);%theta
x(4,1)=Ci(4);
ua(1)=0;

for i=2:1:pasos
    estado=[x(1,i-1);x(2,i-1);x(3,i-1);x(4,i-1)];%Guardo el estado
    integracion=x(4,i-1)+delta_t*(ref(1,i-1)-Mat_Cc*estado);%Es el termino de la pseudo inversa
    u_actual=-Ka*estado;
    ua=[ua u_actual];
    
    %Ecuaciones del motor
    ia_p=(-Ra/Laa)*estado(1)-(Km/Laa)*estado(2)+(1/Laa)*u_actual;
    w_p=(Ki/J)*estado(1)-(B/J)*estado(2);
    theta_p=estado(2);
    
    xp_actual=[ia_p; w_p; theta_p];
    
    xsig=estado(1:3)+delta_t*xp_actual;
    x(1,i)=xsig(1);
    x(2,i)=xsig(2);
    x(3,i)=xsig(3);
    x(4,i)=integracion;
end

figure
plot(t,x(3,1:end));title('Angulo Theta y referencia');hold on;grid on;
plot(t,ref);legend('Angulo Theta','ref');

figure
plot(t,x(2,:));title('Velocidad angular, w');hold on;

figure
plot(t,x(1,:));title('Corriente Ia');hold on;

figure
plot(t,ua);title('Accion de control');hold on;

figure
plot(t,T11);title('Torque');hold on;