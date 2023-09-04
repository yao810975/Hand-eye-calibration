clear 
clc
% Noise-free test code
%% parameter setting
% Set the hand-eye relationship Tc_b ，p_e
Tc_e =  Td([56,17,15])*R(1,0.5,1.1);
p_b = [30, 32 ,510];

% UR5机器人参数 
DH = [  89.2     0        pi/2
        0       -425        0
        0       -392        0
        109.3    0        pi/2
        9.475    0       -pi/2
        8.25     0           0];
   
    
%% Simulation to obtain unbiased measurement data
% p1, p2, and p3 are set for the rotation of joint 6, used to determine the fitting axis
% pc is not coplanar with p1, p2, and p3, and is used for the final solution
q1 = [2   -129.16 -25.67 -85.39 0 66.66];
q2 = [2  -129.16 -25.67 -85.39 0 96.66];
q3 = [2  -129.16 -25.67 -85.39 0 266.66];
qc = [0   -129.16 -25.67 -85.39 0 86.66];

% Tb_e
[T1] = fkine(DH,q1*pi/180);      
[T2] = fkine(DH,q2*pi/180);      
[T3] = fkine(DH,q3*pi/180); 
[Tc] = fkine(DH,qc*pi/180);

% p_c = Tc_e*Te_b*p_b 
p1 = Tc_e*inv(T1)*[p_b,1]';p2 = Tc_e*inv(T2)*[p_b,1]';
p3 = Tc_e*inv(T3)*[p_b,1]';pc = Tc_e*inv(Tc)*[p_b,1]';
p1 = [p1(1) p1(2) p1(3)]; p2 = [p2(1) p2(2) p2(3)]; 
p3 = [p3(1) p3(2) p3(3)]; pc = [pc(1) pc(2) pc(3)];
P_cir = [p1;p2;p3];

% Take any point in p1, p2, p3 and pc to form the solution point set P_cal
P_cal  = [p1;pc];
T(:,:,1) = T1;
T(:,:,2) = Tc;

%% Start solving: Step1
% Tc_e=Dc_e(Tdx1,Tdy,Tdz)*Rc_e(thetax1,thetay,thetaz);
[f1,c1,r1] = fit_circle(P_cir);
% If the error is too large, consider the negative direction f=-f
[thetax,thetay] = thetaxy(f1');

%% Solution 1:

k = [];
for i = 1:size(P_cal,1)
    [ki,mi,wi,kyi] = coefficient_1(inv(T(:,:,i)),thetax,thetay,P_cal(i,:),f1,c1);
    K(:,:,3*i-2) = ki;
    K(:,:,3*i-1) = mi;
    K(:,:,3*i) = wi;
    k = [k;kyi];
end
% x(1) = Tdy; x(2)=p_e_x;   x(3)=p_e_z;  x(4)=p_e_y; x(5) = thetaz ;
x = fsolve(@(x) function_1(x,K,k,3*size(P_cal,1)),[0;0;0;0;0],optimoptions('fsolve','Display','iter'));
p_b_1 = [x(2),x(4),x(3)];

% 根据直线方程求解Tdx，Tdz
Tdy = x(1);
Tdx = f1(1)*(Tdy-c1(2))/f1(2)+c1(1);
Tdz = f1(3)*(Tdy-c1(2))/f1(2)+c1(3);
Tc_e_1 = Td([Tdx,Tdy,Tdz])*R(thetax,thetay,x(5));
%% SoLution 2:

for i = 1:size(P_cal,1)
    [kxi,kyi,kzi] = coefficient_2(T(:,:,i),thetax,thetay,P_cal(i,:),f1,c1);
    kx(:,:,i) = kxi;
    ky(:,:,i) = kyi;
    kz(:,:,i) = kzi;
end

%  x2(1) = thetaz ; x2(2) = Tdy; 
x2 = fsolve(@(x2) function_2(x2,kx,ky,kz,size(P_cal,1)),[0;0],optimoptions('fsolve','Display','iter'));
Tdy = x2(2);
Tdx = f1(1)*(Tdy-c1(2))/f1(2)+c1(1);
Tdz = f1(3)*(Tdy-c1(2))/f1(2)+c1(3);
Tc_b_2 = Td([Tdx,Tdy,Tdz])*R(thetax,thetay,x2(1));

%% verification
p_b_1 - p_b
Tc_b_2 - Tc_e
Tc_e_1 - Tc_e








