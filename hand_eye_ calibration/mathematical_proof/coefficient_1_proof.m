clear
clc
% % % % % % coefficient correctness verification code
%% Parameter setting 
%% known parameter

syms t11 t12 t13 t14 t21 t22 t23 t24 t31 t32 t33 t34 p1 p2 p3 real
%p_c is the measured value, Tb_e can be calculated by the forward kinematics 
% or read by the teaching device.
Tb_e = [t11 t12 t13 t14
       t21 t22 t23 t24
       t31 t32 t33 t34
       0    0   0   1 ]; 
p_c = [p1, p2, p3]; 

% f_c,c_c,thetax,thetay,obtained by Step 1.
syms  f1 f2 f3 c1 c2 c3 thetax thetay real
f_c = [f1, f2,f3];c_c = [c1,c2,c3];

%% unknown parameter£ºy1,thetaz,x2,y2,z2
syms x1 y1 z1 thetaz x2 y2 z2 real

% It can be obtained according to formula (5) in the paper.
x1 = f1*(y1-c2)/f2+c1;
z1 = f3*(y1-c2)/f2+c3;
%  Hand-eye relation matrix parameters, Dmov is translation, Rrot is rotation
Tc_b = Dmov([x1,y1,z1])*Rrot(thetax,thetay,thetaz);
p_e = [x2,y2,z2];


% The calculated value of p_c can be obtained from formula (7)
p_c_x = Tc_b*Tb_e*[p_e 1]';

%% Coefficient of separation:p_c_x
% x direction:
% px = (k1*x + k2*y + k3*z + k4)*cos(thetaz) + ( k5*x + k6*y + k7*z + k8)*sin(thetaz) + k9*x + k10*y + k11*z +  k13*Tdy + k12
k1 = t11*cos(thetay);  k3 = t13*cos(thetay);    k2 = t12*cos(thetay);  k4 = t14*cos(thetay);
k5 = -t21*cos(thetay); k7 = -t23*cos(thetay);   k6 = -t22*cos(thetay); k8 = -t24*cos(thetay);
k9 = t31*sin(thetay);  k11 = t33*sin(thetay);   k10 = t32*sin(thetay); k12 = t34*sin(thetay) + (c1*f2 - c2*f1)/f2;
k13 = f1/f2;

% y direction:
m1 = t21*cos(thetax) + t11*sin(thetax)*sin(thetay); m3 = t23*cos(thetax) + t13*sin(thetax)*sin(thetay); m2 = t22*cos(thetax) + t12*sin(thetax)*sin(thetay);
m4 = t24*cos(thetax) + t14*sin(thetax)*sin(thetay); m5 = t11*cos(thetax) - t21*sin(thetax)*sin(thetay); m7 = t13*cos(thetax) - t23*sin(thetax)*sin(thetay);
m6 = t12*cos(thetax) - t22*sin(thetax)*sin(thetay); m8 = t14*cos(thetax) - t24*sin(thetax)*sin(thetay); m9 = -t31*cos(thetay)*sin(thetax);
m11 = -t33*cos(thetay)*sin(thetax);                 m10 = -t32*cos(thetay)*sin(thetax);                 m12 = -t34*cos(thetay)*sin(thetax);
m13 = 1;

% z direction:
w1 = t21*sin(thetax) - t11*cos(thetax)*sin(thetay); w3 = t23*sin(thetax) - t13*cos(thetax)*sin(thetay); w2 = t22*sin(thetax) - t12*cos(thetax)*sin(thetay);
w4 = t24*sin(thetax) - t14*cos(thetax)*sin(thetay); w5 = t11*sin(thetax) + t21*cos(thetax)*sin(thetay); w7 = t13*sin(thetax) + t23*cos(thetax)*sin(thetay);
w6 = t12*sin(thetax) + t22*cos(thetax)*sin(thetay); w8 = t14*sin(thetax) + t24*cos(thetax)*sin(thetay); w9 = t31*cos(thetax)*cos(thetay);
w11 = t33*cos(thetax)*cos(thetay);                  w10 = t32*cos(thetax)*cos(thetay);                  w12 = t34*cos(thetax)*cos(thetay) - (c2*f3 - c3*f2)/f2;
w13 = f3/f2;

K = [k1 k2 k3 k4;
     k5 k6 k7 k8;
     k9 k10 k11 k12+k13*y1;];
 
M = [m1 m2 m3 m4;
     m5 m6 m7 m8;
     m9 m10 m11 m12+m13*y1;];
 
W = [w1 w2 w3 w4;
     w5 w6 w7 w8;
     w9 w10 w11 w12+w13*y1;];
 %% Verify the correctness of coefficient separation (unknown thetaz x2 y2 z2 y1)
simplify([cos(thetaz) sin(thetaz) 1]*K*[x2; y2; z2; 1]-p_c_x(1))
simplify([cos(thetaz) sin(thetaz) 1]*M*[x2; y2; z2; 1]-p_c_x(2))
simplify([cos(thetaz) sin(thetaz) 1]*W*[x2; y2; z2; 1]-p_c_x(3))


 %% Basic rotation translation function
function [R] = Rrot(thetax,thetay,thetaz)
%Rotation transformation
Rx = [    1           0               0   0
          0    cos(thetax)  -sin(thetax)  0
          0    sin(thetax)   cos(thetax)  0
          0           0               0   1];
      
Ry = [cos(thetay)     0      sin(thetay)  0
          0           1               0   0
      -sin(thetay)    0      cos(thetay)  0
          0           0               0   1];
      
Rz = [cos(thetaz)    -sin(thetaz)     0   0
      sin(thetaz)    cos(thetaz)      0   0
          0         0               1   0
          0         0               0   1];
R = Rx*Ry*Rz;

end

function [T] = Dmov(x)
% translation transformation 


T  = [1         0         0        x(1)
      0         1         0        x(2)
      0         0         1        x(3)
      0         0         0        1];
  

end

