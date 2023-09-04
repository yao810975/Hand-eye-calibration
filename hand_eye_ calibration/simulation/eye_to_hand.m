%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eye-to-hand calibration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% Pcir      :  - Measurement 3D data (fit the circle)
% Pcal      :  - Measurement 3D data
% qc        :  - The robot joint Angle corresponding to Pcal is used to calculate robot pose TB_E
% DH        :  - DH parameters of robot (To calculate robot pose Tb_e)

% OUTPUT
% Tc_b_1        :  - Result of calibration by solution 1
% Tc_b_2        :  - Result of calibration by solution 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Tc_b_1,Tc_b_2]=eye_to_hand(DH,P,P1,q)
% Calculate the robot pose       
for i = 1:size(q,1)
    [Ti] = fkine(DH,q(i,:));
    T(:,:,i) = Ti;
end    
% fit Circle 
[f1,c1,~] = fit_circle(P); 
f1 = -f1';
c1 = c1';
% Solve the rotation parameters thetax, thetay
[thetax,thetay] = thetaxy(f1);

%% Solution 1:Solve Tc_b and p_b
% Gets the coefficient matrix
k = [];
for i = 1:size(P1,1)
    [ki,mi,wi,kyi] = coefficient_1(T(:,:,i),thetax,thetay,P1(i,:),f1,c1);
    K(:,:,3*i-2) = ki;
    K(:,:,3*i-1) = mi;
    K(:,:,3*i) = wi;
    k = [k;kyi];
    
end

% x1(1) = Tdy; x1(2)=p_e_x1;   x1(3)=p_e_z;  x1(4)=p_e_y; x1(5) = thetaz ;
x1 = fsolve(@(x1) function_1(x1,K,k,3*size(P1,1)),[0;0;0;0;0],optimoptions('fsolve','Display','iter'));

p_e_1 = [x1(2),x1(4),x1(3)];
Te_tool = Td(p_e_1);
% Solve x1 and z1 according to the linear equation
Tdy = x1(1);
Tdx1 = f1(1)*(Tdy-c1(2))/f1(2)+c1(1);
Tdz = f1(3)*(Tdy-c1(2))/f1(2)+c1(3);

Tc_b_1 = Td([Tdx1,Tdy,Tdz])*R(thetax,thetay,x1(5));

%% Solution 2:Solve Tc_b only
% Gets the coefficient matrix
for i = 1:size(P1,1)
    [kxi,kyi,kzi] = coefficient_2(inv(T(:,:,i)),thetax,thetay,P1(i,:),f1,c1);
    kx(:,:,i) = kxi;
    ky(:,:,i) = kyi;
    kz(:,:,i) = kzi;
end

%  x2(1) = thetaz ; x2(2) = y1; 
x2 = fsolve(@(x2) function_2(x2,kx,ky,kz,size(P1,1)),[0;0],optimoptions('fsolve','Display','iter'));
Tdy = x2(2);
Tdx2 = f1(1)*(Tdy-c1(2))/f1(2)+c1(1);
Tdz = f1(3)*(Tdy-c1(2))/f1(2)+c1(3);
Tc_b_2 = Td([Tdx2,Tdy,Tdz])*R(thetax,thetay,x2(1));


% save('Tc_b_1.txt','Tc_b_1','-ascii')
% save('Tc_b_2.txt','Tc_b_2','-ascii')

end


