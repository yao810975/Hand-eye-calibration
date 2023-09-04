clear 
clc
% Noise-free test code(The method requires at least four
% points to complete the solution)
%% parameter setting
% Set the hand-eye relationship Tc_b £¬p_e
Tc_b = Td([300,-200,560])*R(0.421,0.120,-0.121);
p_e = [10, 6.23 ,9.12];

% UR5 robot parameters
DH = [  89.2     0        pi/2
        0       -425        0
        0       -392        0
        109.3    0        pi/2
        9.475    0       -pi/2
        8.25     0           0];
   
    
%% Simulation to obtain unbiased measurement data
% p1, p2, and p3 are set for the rotation of joint 1, used to determine the fitting axis
% pc is not coplanar with p1, p2, and p3, and is used for the final solution
q1 = [4   -129.16 -25.67 -85.39 0 86.66];
q2 = [2   -129.16 -25.67 -85.39 0 86.66];
q3 = [1   -129.16 -25.67 -85.39 0 86.66];

qc = [0   -10      -05.67 -2   0.1 16.66];

% Tb_e
[T1] = fkine(DH,q1*pi/180);      
[T2] = fkine(DH,q2*pi/180);      
[T3] = fkine(DH,q3*pi/180); 
[Tc] = fkine(DH,qc*pi/180);

% p_c =Tc_b*Tb_e*P_e
p1 =Tc_b*T1*[p_e,1]'; p2 = Tc_b*T2*[p_e,1]';
p3 = Tc_b*T3*[p_e,1]';pc = Tc_b*Tc*[p_e,1]';
p1 = [p1(1) p1(2) p1(3)]; p2 = [p2(1) p2(2) p2(3)]; 
p3 = [p3(1) p3(2) p3(3)]; pc = [pc(1) pc(2) pc(3)];
P_cir = [p1;p2;p3];

% Take any point in p1, p2, p3 and pc to form the solution point set P_cal
P_cal  = [p1;pc];

T(:,:,1) = T1;
T(:,:,2) = Tc;

%% Start solving: Step1
% Tc_b=Dc_b(Tdx1,Tdy,Tdz)*Rc_b(thetax1,thetay,thetaz);
% Fitting circle
[f1,c1,r1] = fit_circle(P_cir);
% If the error is too large, consider the negative direction f=-f
[thetax,thetay] = thetaxy(f1');

%% Solution 1:

k = [];
for i = 1:size(P_cal,1)
    [ki,mi,wi,kyi] = coefficient_1(T(:,:,i),thetax,thetay,P_cal(i,:),f1,c1);
    K(:,:,3*i-2) = ki;
    K(:,:,3*i-1) = mi;
    K(:,:,3*i) = wi;
    k = [k;kyi];
    
end

% x1(1) = ; x1(2)=p_e_x1;   x1(3)=p_e_z;  x1(4)=p_e_y; x1(5) = thetaz ;
x1 = fsolve(@(x1) function_1(x1,K,k,3*size(P_cal,1)),[0;0;0;0;0],optimoptions('fsolve','Display','iter'));

p_e_1 = [x1(2),x1(4),x1(3)];
% Solve x1, z1 according to the linear equation
Tdy = x1(1);
Tdx = f1(1)*(Tdy-c1(2))/f1(2)+c1(1);
Tdz = f1(3)*(Tdy-c1(2))/f1(2)+c1(3);
Tc_b_1 = Td([Tdx,Tdy,Tdz])*R(thetax,thetay,x1(5));

%% SoLution 2:
for i = 1:size(P_cal,1)
    [kxi,kyi,kzi] = coefficient_2(inv(T(:,:,i)),thetax,thetay,P_cal(i,:),f1,c1);
    kx(:,:,i) = kxi;
    ky(:,:,i) = kyi;
    kz(:,:,i) = kzi;
end

%  x2(1) = thetaz ; x2(2) = y1; 
x2 = fsolve(@(x2) function_2(x2,kx,ky,kz,size(P_cal,1)),[0;0],optimoptions('fsolve','Display','iter'));
Tdy = x2(2);
Tdx = f1(1)*(Tdy-c1(2))/f1(2)+c1(1);
Tdz = f1(3)*(Tdy-c1(2))/f1(2)+c1(3);
Tc_b_2 = Td([Tdx,Tdy,Tdz])*R(thetax,thetay,x2(1));

%% verification
Tc_b - Tc_b_1
p_e_1- p_e
Tc_b - Tc_b_2




