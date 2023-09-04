%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eye-in-hand calibration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% Pcir      :  - Measurement 3D data (fit the circle)
% Pcal      :  - Measurement 3D data
% qc        :  - The robot joint Angle corresponding to Pcal is used to calculate robot pose TB_E
% DH        :  - DH parameters of robot (To calculate robot pose Tb_e)

% OUTPUT
% Tc_e_1        :  - Result of calibration by solution 1
% Tc_e_2        :  - Result of calibration by solution 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Tc_e_1,Tc_e_2,f1,c1]=eye_in_hand(DH,Pcir,Pcal,qc)
 
 % Calculate the robot pose      
for i = 1:size(qc,1)
    Ti = fkine(DH,qc(i,:));
    T_2(:,:,i) = Ti;
    T_1(:,:,i) = inv(Ti);
end

% fit Circle 
[f1,c1,~] = fit_circle(Pcir);
% Solve the rotation parameters thetax, thetay
[thetax,thetay] = thetaxy(f1');

%% Solution 1:Solve Tc_e and p_b
% Gets the coefficient matrix
k = [];
for i = 1:size(Pcal,1)
    [ki,mi,wi,kyi] = coefficient_1(T_1(:,:,i),thetax,thetay,Pcal(i,:),f1,c1);
    K(:,:,3*i-2) = ki;
    K(:,:,3*i-1) = mi;
    K(:,:,3*i) = wi;
    k = [k;kyi];
    
end

% x(1) = Tdy; x(2)=p_e_x;   x(3)=p_e_z;  x(4)=p_e_y; x(5) = thetaz ;
% options = optimoptions('fsolve','Display','iter');
options = optimoptions('fsolve');
x1 = fsolve(@(x) function_1(x,K,k,3*size(Pcal)),[0;0;0;0;0],options);
p_b_1 = [x1(2),x1(4),x1(3)];

% Solve Tdx and Tdz according to the linear equation
Tdy_1 = x1(1);
Tdx_1 = f1(1)*(Tdy_1-c1(2))/f1(2)+c1(1);
Tdz_1 = f1(3)*(Tdy_1-c1(2))/f1(2)+c1(3);
Tc_e_1 = Td([Tdx_1,Tdy_1,Tdz_1])*R(thetax,thetay,x1(5));


%% Solution 2:Solve Tc_e only
% 
for i = 1:size(Pcal,1)
    [kxi,kyi,kzi] = coefficient_2(T_2(:,:,i),thetax,thetay,Pcal(i,:),f1,c1);
    kx(:,:,i) = kxi;
    ky(:,:,i) = kyi;
    kz(:,:,i) = kzi;
end


x2 = fsolve(@(x) function_2(x,kx,ky,kz,size(Pcal,1)),[0;0],options);

Tdy_2 = x2(2);
Tdx_2 = f1(1)*(Tdy_2-c1(2))/f1(2)+c1(1);
Tdz_2 = f1(3)*(Tdy_2-c1(2))/f1(2)+c1(3);

Tc_e_2 = R(thetax,thetay,x2(1));
Tc_e_2(1,4) = Tdx_2;
Tc_e_2(2,4) = Tdy_2;
Tc_e_2(3,4) = Tdz_2;
% save('.\Tc_e_1.txt','Tc_e_1','-ascii')
% save('.\Tc_e_2.txt','Tc_e_2','-ascii')
end



