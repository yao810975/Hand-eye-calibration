clear
clc
% % % % ��ⷽʽ��ϵ����ȷ����֤����
%% �������ã�δ����Ϊthetaz x2 y2 z2 y1)
syms t11 t12 t13 t14 t21 t22 t23 t24 t31 t32 t33 t34  real
% Tt_b���������Tool��Baseת������
Tt_b = [t11 t12 t13 t14
       t21 t22 t23 t24
       t31 t32 t33 t34
       0    0   0   1 ]; 

syms  f1 f2 f3 c1 c2 c3 p1 p2 p3 x2 y2 z2 real
% f_camΪ���Բƽ�淨����������ȡ�ؽ�1z�������򣩣�c_camΪ���Բ��Բ��,p_tΪ��������ϵ�±�־��λ��
f_cam = [f1, f2,f3];
c_cam = [c1,c2,c3];
p_c = [p1, p2, p3];
p_t = [x2,y2,z2];

% �������۹�ϵ����ɿ�����ƽ������ת{Base}��{Cam}�ɱ�ʾΪ
syms x1 y1 z1 thetax thetay thetaz
Tc_b = Dmov([x1,y1,z1])*Rrot(thetax,thetay,thetaz);
% �������й�ʽ(5)
x1 = f1*(y1-c2)/f2+c1;
z1 = f3*(y1-c2)/f2+c3;
% �ٴδ���
Tc_b = Dmov([x1,y1,z1])*Rrot(thetax,thetay,thetaz);

% �����Ĺ�ʽ(13),�ɵ�p_t����ֵ
p_t_x = Tt_b*inv(Tc_b)*[p_c 1]';
%% ����p_t_x��ϵ��

% ��p_t_x����ϵ����δ֪��Ϊthetaz,y1)
% (k1*y1 + k2)*cos(thetaz) + (k3*y1 + k4)*sin(thetaz) + k5*y1 + k6 - p_t_x(1)*f2
k1 = f3*t11*cos(thetax)*sin(thetay) - f2*t12*cos(thetax) - f3*t12*sin(thetax) - f2*t11*sin(thetax)*sin(thetay) - f1*t11*cos(thetay);
k2 = c2*f1*t11*cos(thetay) - c1*f2*t11*cos(thetay) + f2*p1*t11*cos(thetay) + f2*p2*t12*cos(thetax) + c2*f3*t12*sin(thetax) - c3*f2*t12*sin(thetax) + f2*p3*t12*sin(thetax) - c2*f3*t11*cos(thetax)*sin(thetay) + c3*f2*t11*cos(thetax)*sin(thetay) - f2*p3*t11*cos(thetax)*sin(thetay) + f2*p2*t11*sin(thetax)*sin(thetay);
k3 = f1*t12*cos(thetay) - f2*t11*cos(thetax) - f3*t11*sin(thetax) + f2*t12*sin(thetax)*sin(thetay) - f3*t12*cos(thetax)*sin(thetay);
k4 = c1*f2*t12*cos(thetay) - c2*f1*t12*cos(thetay) - f2*p1*t12*cos(thetay) + f2*p2*t11*cos(thetax) + c2*f3*t11*sin(thetax) - c3*f2*t11*sin(thetax) + f2*p3*t11*sin(thetax) + c2*f3*t12*cos(thetax)*sin(thetay) - c3*f2*t12*cos(thetax)*sin(thetay) + f2*p3*t12*cos(thetax)*sin(thetay) - f2*p2*t12*sin(thetax)*sin(thetay);
k5 = f2*t13*cos(thetay)*sin(thetax) - f3*t13*cos(thetax)*cos(thetay) - f1*t13*sin(thetay);
k6 = f2*t14 - c1*f2*t13*sin(thetay) + c2*f1*t13*sin(thetay) + f2*p1*t13*sin(thetay) + c2*f3*t13*cos(thetax)*cos(thetay) - c3*f2*t13*cos(thetax)*cos(thetay) + f2*p3*t13*cos(thetax)*cos(thetay) - f2*p2*t13*cos(thetay)*sin(thetax);

m1 = f3*t21*cos(thetax)*sin(thetay) - f2*t22*cos(thetax) - f3*t22*sin(thetax) - f2*t21*sin(thetax)*sin(thetay) - f1*t21*cos(thetay);
m2 = c2*f1*t21*cos(thetay) - c1*f2*t21*cos(thetay) + f2*p1*t21*cos(thetay) + f2*p2*t22*cos(thetax) + c2*f3*t22*sin(thetax) - c3*f2*t22*sin(thetax) + f2*p3*t22*sin(thetax) - c2*f3*t21*cos(thetax)*sin(thetay) + c3*f2*t21*cos(thetax)*sin(thetay) - f2*p3*t21*cos(thetax)*sin(thetay) + f2*p2*t21*sin(thetax)*sin(thetay);
m3 = f1*t22*cos(thetay) - f2*t21*cos(thetax) - f3*t21*sin(thetax) + f2*t22*sin(thetax)*sin(thetay) - f3*t22*cos(thetax)*sin(thetay);
m4 = c1*f2*t22*cos(thetay) - c2*f1*t22*cos(thetay) - f2*p1*t22*cos(thetay) + f2*p2*t21*cos(thetax) + c2*f3*t21*sin(thetax) - c3*f2*t21*sin(thetax) + f2*p3*t21*sin(thetax) + c2*f3*t22*cos(thetax)*sin(thetay) - c3*f2*t22*cos(thetax)*sin(thetay) + f2*p3*t22*cos(thetax)*sin(thetay) - f2*p2*t22*sin(thetax)*sin(thetay);
m5 = f2*t23*cos(thetay)*sin(thetax) - f3*t23*cos(thetax)*cos(thetay) - f1*t23*sin(thetay);
m6 = f2*t24 - c1*f2*t23*sin(thetay) + c2*f1*t23*sin(thetay) + f2*p1*t23*sin(thetay) + c2*f3*t23*cos(thetax)*cos(thetay) - c3*f2*t23*cos(thetax)*cos(thetay) + f2*p3*t23*cos(thetax)*cos(thetay) - f2*p2*t23*cos(thetay)*sin(thetax);


w1 = f3*t31*cos(thetax)*sin(thetay) - f2*t32*cos(thetax) - f3*t32*sin(thetax) - f2*t31*sin(thetax)*sin(thetay) - f1*t31*cos(thetay);
w2 = c2*f1*t31*cos(thetay) - c1*f2*t31*cos(thetay) + f2*p1*t31*cos(thetay) + f2*p2*t32*cos(thetax) + c2*f3*t32*sin(thetax) - c3*f2*t32*sin(thetax) + f2*p3*t32*sin(thetax) - c2*f3*t31*cos(thetax)*sin(thetay) + c3*f2*t31*cos(thetax)*sin(thetay) - f2*p3*t31*cos(thetax)*sin(thetay) + f2*p2*t31*sin(thetax)*sin(thetay);
w3 = f1*t32*cos(thetay) - f2*t31*cos(thetax) - f3*t31*sin(thetax) + f2*t32*sin(thetax)*sin(thetay) - f3*t32*cos(thetax)*sin(thetay);
w4 = c1*f2*t32*cos(thetay) - c2*f1*t32*cos(thetay) - f2*p1*t32*cos(thetay) + f2*p2*t31*cos(thetax) + c2*f3*t31*sin(thetax) - c3*f2*t31*sin(thetax) + f2*p3*t31*sin(thetax) + c2*f3*t32*cos(thetax)*sin(thetay) - c3*f2*t32*cos(thetax)*sin(thetay) + f2*p3*t32*cos(thetax)*sin(thetay) - f2*p2*t32*sin(thetax)*sin(thetay);
w5 = f2*t33*cos(thetay)*sin(thetax) - f3*t33*cos(thetax)*cos(thetay) - f1*t33*sin(thetay);
w6 = f2*t34 - c1*f2*t33*sin(thetay) + c2*f1*t33*sin(thetay) + f2*p1*t33*sin(thetay) + c2*f3*t33*cos(thetax)*cos(thetay) - c3*f2*t33*cos(thetax)*cos(thetay) + f2*p3*t33*cos(thetax)*cos(thetay) - f2*p2*t33*cos(thetay)*sin(thetax);

K = [k1 k2 
     k3 k4 
     k5 k6 ];
K = K/f2 ;

M = [m1 m2 
     m3 m4
     m5 m6 ];
M = M/f2;

W = [w1 w2 
     w3 w4
     w5 w6];
W = W/f2;

 %% ��֤ϵ��������ȷ�ԣ�δ֪��thetaz  y1)
simplify([cos(thetaz) sin(thetaz) 1]*K*[y1;1]-p_t_x(1))
simplify([cos(thetaz) sin(thetaz) 1]*M*[y1;1]-p_t_x(2))
simplify([cos(thetaz) sin(thetaz) 1]*W*[y1;1]-p_t_x(3))


 %% ������תƽ�ƺ���
function [R] = Rrot(thetax,thetay,thetaz)
%��ת�任����
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
%ƽ�Ʊ任����


T  = [1         0         0        x(1)
      0         1         0        x(2)
      0         0         1        x(3)
      0         0         0        1];
  

end

