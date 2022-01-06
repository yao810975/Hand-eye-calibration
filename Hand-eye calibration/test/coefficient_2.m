%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p_e = Tb_e*Te_c*p_c
% p_e = Te_b*Tb_c*p_c
% 分离常系数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% T        :  - eye_in_hand输入Tb_e;eye_to_hand输入Te_b
% thetax   :  - 绕x旋转的弧度值
% thetay   :  - 绕y旋转的弧度值
% p        :  - 相机测量点
% f        :  - 圆平面法向量
% c        :  - 圆心
% OUTPUT
% k        :  - x方向常系数
% m        :  - y方向常系数
% w        :  - z方向常系数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [k,m,w] = coefficient_2(T,thetax,thetay,p,f,c)

t11 = T(1,1);t12 = T(1,2);t13 = T(1,3);t14 = T(1,4);
t21 = T(2,1);t22 = T(2,2);t23 = T(2,3);t24 = T(2,4);
t31 = T(3,1);t32 = T(3,2);t33 = T(3,3);t34 = T(3,4);
f1x = f(1); f1y = f(2); f1z = f(3);
c1x = c(1); c1y = c(2); c1z = c(3);
pcx = p(1); pcy = p(2); pcz = p(3);


%  p_e_x = nx/f1y 
%  dx = f1y
k1 = f1z*t11*cos(thetax)*sin(thetay) - f1y*t12*cos(thetax) - f1z*t12*sin(thetax) - f1y*t11*sin(thetax)*sin(thetay) - f1x*t11*cos(thetay);
k2 = c1y*f1x*t11*cos(thetay) - c1x*f1y*t11*cos(thetay) + f1y*pcx*t11*cos(thetay) + f1y*pcy*t12*cos(thetax) + c1y*f1z*t12*sin(thetax) - c1z*f1y*t12*sin(thetax) + f1y*pcz*t12*sin(thetax) - c1y*f1z*t11*cos(thetax)*sin(thetay) + c1z*f1y*t11*cos(thetax)*sin(thetay) - f1y*pcz*t11*cos(thetax)*sin(thetay) + f1y*pcy*t11*sin(thetax)*sin(thetay);
k3 = f1x*t12*cos(thetay) - f1y*t11*cos(thetax) - f1z*t11*sin(thetax) + f1y*t12*sin(thetax)*sin(thetay) - f1z*t12*cos(thetax)*sin(thetay);
k4 = c1x*f1y*t12*cos(thetay) - c1y*f1x*t12*cos(thetay) - f1y*pcx*t12*cos(thetay) + f1y*pcy*t11*cos(thetax) + c1y*f1z*t11*sin(thetax) - c1z*f1y*t11*sin(thetax) + f1y*pcz*t11*sin(thetax) + c1y*f1z*t12*cos(thetax)*sin(thetay) - c1z*f1y*t12*cos(thetax)*sin(thetay) + f1y*pcz*t12*cos(thetax)*sin(thetay) - f1y*pcy*t12*sin(thetax)*sin(thetay);
k5 = f1y*t13*cos(thetay)*sin(thetax) - f1z*t13*cos(thetax)*cos(thetay) - f1x*t13*sin(thetay);
k6 = f1y*t14 - c1x*f1y*t13*sin(thetay) + c1y*f1x*t13*sin(thetay) + f1y*pcx*t13*sin(thetay) + c1y*f1z*t13*cos(thetax)*cos(thetay) - c1z*f1y*t13*cos(thetax)*cos(thetay) + f1y*pcz*t13*cos(thetax)*cos(thetay) - f1y*pcy*t13*cos(thetay)*sin(thetax);

m1 = f1z*t21*cos(thetax)*sin(thetay) - f1y*t22*cos(thetax) - f1z*t22*sin(thetax) - f1y*t21*sin(thetax)*sin(thetay) - f1x*t21*cos(thetay);
m2 = c1y*f1x*t21*cos(thetay) - c1x*f1y*t21*cos(thetay) + f1y*pcx*t21*cos(thetay) + f1y*pcy*t22*cos(thetax) + c1y*f1z*t22*sin(thetax) - c1z*f1y*t22*sin(thetax) + f1y*pcz*t22*sin(thetax) - c1y*f1z*t21*cos(thetax)*sin(thetay) + c1z*f1y*t21*cos(thetax)*sin(thetay) - f1y*pcz*t21*cos(thetax)*sin(thetay) + f1y*pcy*t21*sin(thetax)*sin(thetay);
m3 = f1x*t22*cos(thetay) - f1y*t21*cos(thetax) - f1z*t21*sin(thetax) + f1y*t22*sin(thetax)*sin(thetay) - f1z*t22*cos(thetax)*sin(thetay);
m4 = c1x*f1y*t22*cos(thetay) - c1y*f1x*t22*cos(thetay) - f1y*pcx*t22*cos(thetay) + f1y*pcy*t21*cos(thetax) + c1y*f1z*t21*sin(thetax) - c1z*f1y*t21*sin(thetax) + f1y*pcz*t21*sin(thetax) + c1y*f1z*t22*cos(thetax)*sin(thetay) - c1z*f1y*t22*cos(thetax)*sin(thetay) + f1y*pcz*t22*cos(thetax)*sin(thetay) - f1y*pcy*t22*sin(thetax)*sin(thetay);
m5 = f1y*t23*cos(thetay)*sin(thetax) - f1z*t23*cos(thetax)*cos(thetay) - f1x*t23*sin(thetay);
m6 = f1y*t24 - c1x*f1y*t23*sin(thetay) + c1y*f1x*t23*sin(thetay) + f1y*pcx*t23*sin(thetay) + c1y*f1z*t23*cos(thetax)*cos(thetay) - c1z*f1y*t23*cos(thetax)*cos(thetay) + f1y*pcz*t23*cos(thetax)*cos(thetay) - f1y*pcy*t23*cos(thetay)*sin(thetax);


w1 = f1z*t31*cos(thetax)*sin(thetay) - f1y*t32*cos(thetax) - f1z*t32*sin(thetax) - f1y*t31*sin(thetax)*sin(thetay) - f1x*t31*cos(thetay);
w2 = c1y*f1x*t31*cos(thetay) - c1x*f1y*t31*cos(thetay) + f1y*pcx*t31*cos(thetay) + f1y*pcy*t32*cos(thetax) + c1y*f1z*t32*sin(thetax) - c1z*f1y*t32*sin(thetax) + f1y*pcz*t32*sin(thetax) - c1y*f1z*t31*cos(thetax)*sin(thetay) + c1z*f1y*t31*cos(thetax)*sin(thetay) - f1y*pcz*t31*cos(thetax)*sin(thetay) + f1y*pcy*t31*sin(thetax)*sin(thetay);
w3 = f1x*t32*cos(thetay) - f1y*t31*cos(thetax) - f1z*t31*sin(thetax) + f1y*t32*sin(thetax)*sin(thetay) - f1z*t32*cos(thetax)*sin(thetay);
w4 = c1x*f1y*t32*cos(thetay) - c1y*f1x*t32*cos(thetay) - f1y*pcx*t32*cos(thetay) + f1y*pcy*t31*cos(thetax) + c1y*f1z*t31*sin(thetax) - c1z*f1y*t31*sin(thetax) + f1y*pcz*t31*sin(thetax) + c1y*f1z*t32*cos(thetax)*sin(thetay) - c1z*f1y*t32*cos(thetax)*sin(thetay) + f1y*pcz*t32*cos(thetax)*sin(thetay) - f1y*pcy*t32*sin(thetax)*sin(thetay);
w5 = f1y*t33*cos(thetay)*sin(thetax) - f1z*t33*cos(thetax)*cos(thetay) - f1x*t33*sin(thetay);
w6 = f1y*t34 - c1x*f1y*t33*sin(thetay) + c1y*f1x*t33*sin(thetay) + f1y*pcx*t33*sin(thetay) + c1y*f1z*t33*cos(thetax)*cos(thetay) - c1z*f1y*t33*cos(thetax)*cos(thetay) + f1y*pcz*t33*cos(thetax)*cos(thetay) - f1y*pcy*t33*cos(thetay)*sin(thetax);

% 返回三个常数K矩阵
k = [k1 k2 
     k3 k4 
     k5 k6 ];
 
m = [m1 m2 
     m3 m4
     m5 m6 ];
 
w = [w1 w2 
     w3 w4
     w5 w6];


end

% 西南科技大学特种机器人实验室2021/9/10
