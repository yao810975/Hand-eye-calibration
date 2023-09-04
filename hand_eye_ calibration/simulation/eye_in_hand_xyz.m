%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eye-in-hand calibration uses only single directional data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% Pcir      :  - Measurement 3D data (fit the circle)
% Pcal      :  - Measurement 3D data
% qc        :  - The robot joint Angle corresponding to Pcal is used to calculate robot pose TB_E
% DH        :  - DH parameters of robot (To calculate robot pose Tb_e)

% OUTPUT
% Tc_e_1_x        :  - Result of calibration by solution 1_x 
% Tc_e_1_y        :  - Result of calibration by solution 1_y
% Tc_e_1_z        :  - Result of calibration by solution 1_z
% Tc_e_2_x        :  - Result of calibration by solution 2_x
% Tc_e_2_y        :  - Result of calibration by solution 2_y
% Tc_e_2_z        :  - Result of calibration by solution 2_z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Tc_e_1_x,Tc_e_1_y,Tc_e_1_z,Tc_e_2_x,Tc_e_2_y,Tc_e_2_z] = eye_in_hand_xyz(DH,Pcir,Pcal,qc)

    for i = 1:size(qc,1)
        Ti = fkine(DH,qc(i,:));
        T_2(:,:,i) = Ti;
        T_1(:,:,i) = inv(Ti);
    end
    
    [f1,c1,~] = fit_circle(Pcir);

    [thetax,thetay] = thetaxy(f1');

    %%  solution 1
    % Gets the coefficient matrix
    k = [];
    for i = 1:length(Pcal)
        [ki,mi,wi,kyi] = coefficient_1(T_1(:,:,i),thetax,thetay,Pcal(i,:),f1,c1);
    % only x direction coefficient
        K_x(:,:,i) = ki;
        k_x(i) = kyi(1);    
    % only y direction coefficient
        K_y(:,:,i) = mi;
        k_y(i) = kyi(2);
    % only z direction coefficient
        K_z(:,:,i) = wi;
        k_z(i) = kyi(3);  
    end
    options = optimoptions('fsolve','Display','iter');
    options = optimoptions(options,'StepTolerance',1e-16);
   %%   x direction
    x1 = fsolve(@(x) function_1(x,K_x,k_x,size(Pcal)),[0;0;0;0;0],options);
    p_b_1 = [x1(2),x1(4),x1(3)];

    Tdy_1 = x1(1);
    Tdx_1 = f1(1)*(Tdy_1-c1(2))/f1(2)+c1(1);
    Tdz_1 = f1(3)*(Tdy_1-c1(2))/f1(2)+c1(3);
    Tc_e_1_x = Td([Tdx_1,Tdy_1,Tdz_1])*R(thetax,thetay,x1(5));
    x1 = fsolve(@(x) function_1(x,K_y,k_y,size(Pcal)),[0;0;0;0;0],options);
    
    p_b_1 = [x1(2),x1(4),x1(3)];
    %%   y direction

    Tdy_1 = x1(1);
    Tdx_1 = f1(1)*(Tdy_1-c1(2))/f1(2)+c1(1);
    Tdz_1 = f1(3)*(Tdy_1-c1(2))/f1(2)+c1(3);
    Tc_e_1_y = Td([Tdx_1,Tdy_1,Tdz_1])*R(thetax,thetay,x1(5));
     %%   z direction
    x1 = fsolve(@(x) function_1(x,K_z,k_z,size(Pcal)),[0;0;0;0;0],options);
    p_b_1 = [x1(2),x1(4),x1(3)];

    Tdy_1 = x1(1);
    Tdx_1 = f1(1)*(Tdy_1-c1(2))/f1(2)+c1(1);
    Tdz_1 = f1(3)*(Tdy_1-c1(2))/f1(2)+c1(3);
    Tc_e_1_z = Td([Tdx_1,Tdy_1,Tdz_1])*R(thetax,thetay,x1(5));

    %% solution 2
    
    for i = 1:length(Pcal)
        [kxi,kyi,kzi] = coefficient_2(T_2(:,:,i),thetax,thetay,Pcal(i,:),f1,c1);
        kx(:,:,i) = kxi;
        ky(:,:,i) = kyi;
        kz(:,:,i) = kzi;
    end
    %% x direction

    x2 = fsolve(@(x) function_2(x,kx,kx,kx,length(Pcal)),[0;0],options);

    Tdy_2 = x2(2);
    Tdx_2 = f1(1)*(Tdy_2-c1(2))/f1(2)+c1(1);
    Tdz_2 = f1(3)*(Tdy_2-c1(2))/f1(2)+c1(3);
    Tc_e_2_x = Td([Tdx_2,Tdy_2,Tdz_2])*R(thetax,thetay,x2(1));
    %% y direction

    x2 = fsolve(@(x) function_2(x,ky,ky,ky,length(Pcal)),[0;0],options);

    Tdy_2 = x2(2);
    Tdx_2 = f1(1)*(Tdy_2-c1(2))/f1(2)+c1(1);
    Tdz_2 = f1(3)*(Tdy_2-c1(2))/f1(2)+c1(3);
    Tc_e_2_y = Td([Tdx_2,Tdy_2,Tdz_2])*R(thetax,thetay,x2(1));
    %% z direction
    x2 = fsolve(@(x) function_2(x,kz,kz,kz,length(Pcal)),[0;0],options);

    Tdy_2 = x2(2);
    Tdx_2 = f1(1)*(Tdy_2-c1(2))/f1(2)+c1(1);
    Tdz_2 = f1(3)*(Tdy_2-c1(2))/f1(2)+c1(3);
    Tc_e_2_z = Td([Tdx_2,Tdy_2,Tdz_2])*R(thetax,thetay,x2(1));


    end