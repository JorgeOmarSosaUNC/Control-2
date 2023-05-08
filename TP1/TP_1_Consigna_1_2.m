% Jorge Omar Sosa
% Sistema de Control II - 2023
% Trabajo Practico 1 - 1.2
% ===========================================================
R=5.6e3; L=10e-6; C=100e-9; vin=12;
T=5.0e-3; tp=0.5e-10; items=T/tp; t=linspace(0, T, items);
Ip=0; Vcp=0;
I=zeros(1, items); Vc=zeros(1,items); u=linspace(0, 0, items);
I(1)=0; Vc(1)=0; u(1)=vin;
MAT_A=[-R/L -1/L ; 1/C 0];
MAT_B=[1/L ; 0];
MAT_C=[R 0];
Il(1)=0; Vcl(1)=0; x=[I(1) Vc(1)]'; y(1)=0; xop=[0 0]';
>> cont=0;
for i=1:items-1
cont = cont + tp;
if (cont>=1e-3)
cont=0;
vin=vin*-1;
end
estado=[I(i); Vc(i)]; u(i)=vin;
xp=MAT_A*(x-xop)+MAT_B*u(i);
x=x+xp*tp;
Y=MAT_C*x;
Il(i+1) = x(1);
Vcl(i+1)= x(2);
i=i+1;
end
figure(1);
subplot(3, 1, 1);
plot(t, Il, 'r'); title('corriente, i_t');
subplot(3, 1, 2);
plot(t, Vcl, 'r'); title('tension Vc, v_t');
subplot(3, 1, 3);
plot(t, u, 'r'); title('tension de entrada, vin_t')