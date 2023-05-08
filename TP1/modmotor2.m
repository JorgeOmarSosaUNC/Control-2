function [X]=modmotor2(t_S, xant, accion, torque) 
Laa=366e-6; 
J=5e-9; 
Ra=55.6; 
B=0; 
Ki=6.49e-3; 
Km=6.53e-3;
TL=torque; 
Va=accion; 
h=1e-7; 
omega=xant(1); 
wp=xant(2); 
ia=xant(3);
for ii=1:t_S/h 
    wpp =(-wp*(Ra*J+Laa*B)-omega*(Ra*B+Ki*Km)+Va*Ki)/(J*Laa); 
    iap=(-Ra*ia-Km*omega+Va)/Laa;
    wp=wp+h*wpp; 
    wp=wp-(TL/J); 
    ia=ia+h*iap;
    omega=omega+h*wp;
end 
X=[omega,wp,ia]; 
end