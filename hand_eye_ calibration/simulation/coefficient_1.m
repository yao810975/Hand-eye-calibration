%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The constant coefficient of solution 1;

% The form of coefficients is     [k1 k2 k3         k4;
%                                  k5 k6 k7         k8;
%                                  k9 k10 k11   k12+k13*y1;];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% T        :  - robot pose: eye_to_hand input Tb_e; eye_in_hand input Te_b;
% thetax   :  - The value of radians about x (solve by 'thetaxy' function)
% thetay   :  - The value of radians about y (solve by 'thetaxy' function)
% p        :  - Measurement 3D data
% f        :  - Normal vector to the circular plane
% c        :  - Center of the circle
% OUTPUT
% k        :  - Constant coefficient in x direction (excluding y1)
% m        :  - Constant coefficient in y direction (excluding y1)
% w        :  - Constant coefficient in z direction (excluding y1)
% ay       :  - Constant coefficient of y1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [k,m,w,ay] = coefficient_1(T,thetax,thetay,p,f,c)

t11 = T(1,1);t12 = T(1,2);t13 = T(1,3);t14 = T(1,4);
t21 = T(2,1);t22 = T(2,2);t23 = T(2,3);t24 = T(2,4);
t31 = T(3,1);t32 = T(3,2);t33 = T(3,3);t34 = T(3,4);
f1x = f(1); f1y = f(2); f1z = f(3);
c1x = c(1); c1y = c(2); c1z = c(3);
% px = (k1*x + k2*z + k3*y + k4)*cos(thetaz) + ( k5*x + k6*z + k7*y + k8)*sin(thetaz) + k9*x + k10*z + k11*y +  k13*Tdy + k12
k1 = t11*cos(thetay);  k2 = t13*cos(thetay);    k3 = t12*cos(thetay);  k4 = t14*cos(thetay);
k5 = -t21*cos(thetay); k6 = -t23*cos(thetay);   k7 = -t22*cos(thetay); k8 = -t24*cos(thetay);
k9 = t31*sin(thetay);  k10 = t33*sin(thetay);   k11 = t32*sin(thetay); k12 = t34*sin(thetay) + (c1x*f1y - c1y*f1x)/f1y - p(1);
k13 = f1x/f1y;

% py = (m1*x + m2*z + m3*y + m4)*cos(thetaz) + ( m5*x + m6*z + m7*y + m8)*sin(thetaz) + m9*x + m10*z + m11*y +  m13*Tdy + m12
m1 = t21*cos(thetax) + t11*sin(thetax)*sin(thetay); m2 = t23*cos(thetax) + t13*sin(thetax)*sin(thetay); m3 = t22*cos(thetax) + t12*sin(thetax)*sin(thetay);
m4 = t24*cos(thetax) + t14*sin(thetax)*sin(thetay); m5 = t11*cos(thetax) - t21*sin(thetax)*sin(thetay); m6 = t13*cos(thetax) - t23*sin(thetax)*sin(thetay);
m7 = t12*cos(thetax) - t22*sin(thetax)*sin(thetay); m8 = t14*cos(thetax) - t24*sin(thetax)*sin(thetay); m9 = -t31*cos(thetay)*sin(thetax);
m10 = -t33*cos(thetay)*sin(thetax);                 m11 = -t32*cos(thetay)*sin(thetax);                 m12 = -t34*cos(thetay)*sin(thetax) - p(2);
m13 = 1;

% pz = (w1*x + w2*z + w3*y + w4)*cos(thetaz) + ( w5*x + w6*z + w7*y + w8)*sin(thetaz) + w9*x + w10*z + w11*y +  w13*Tdy + w12
w1 = t21*sin(thetax) - t11*cos(thetax)*sin(thetay); w2 = t23*sin(thetax) - t13*cos(thetax)*sin(thetay); w3 = t22*sin(thetax) - t12*cos(thetax)*sin(thetay);
w4 = t24*sin(thetax) - t14*cos(thetax)*sin(thetay); w5 = t11*sin(thetax) + t21*cos(thetax)*sin(thetay); w6 = t13*sin(thetax) + t23*cos(thetax)*sin(thetay);
w7 = t12*sin(thetax) + t22*cos(thetax)*sin(thetay); w8 = t14*sin(thetax) + t24*cos(thetax)*sin(thetay); w9 = t31*cos(thetax)*cos(thetay);
w10 = t33*cos(thetax)*cos(thetay);                  w11 = t32*cos(thetax)*cos(thetay);                  w12 = t34*cos(thetax)*cos(thetay) - (c1y*f1z - c1z*f1y)/f1y - p(3);
w13 = f1z/f1y;
% returns three constant K matrices
k = [k1 k2 k3 k4;
     k5 k6 k7 k8;
     k9 k10 k11 k12;];
 
m = [m1 m2 m3 m4;
     m5 m6 m7 m8;
     m9 m10 m11 m12;];
 
w = [w1 w2 w3 w4;
     w5 w6 w7 w8;
     w9 w10 w11 w12;];
 
% ay is y1 constant coefficient (y1 is translation component)
ay = [k13;m13;w13];


end

