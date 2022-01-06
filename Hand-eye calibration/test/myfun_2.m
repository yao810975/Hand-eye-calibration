%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 建立方程组f = p_e - p_e
% p_e = [cos(theta) sin(theta) 1]*K*[Tdy;1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% x        :  - 待求量 x(1) = thetaz ; x(2) = Tdy; 
% kx        :  - x方向系数矩阵p_e_x = [cos(theta) sin(theta) 1]*k*[Tdy;1]
% ky        :  - y方向系数矩阵p_e_y = [cos(theta) sin(theta) 1]*m*[Tdy;1]
% kz        :  - z方向系数矩阵p_e_z = [cos(theta) sin(theta) 1]*w*[Tdy;1]
% K        :  - 由coefficient_2函数结果所组成的系数矩阵
% n        :  - P点的个数

% OUTPUT
% F        :  - 待用fsolve求解的方程组
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = myfun_2(x,kx,ky,kz,n)

for i = 1:n
    p_e_x(i,:) = ([cos(x(1)) sin(x(1)) 1]*kx(:,:,i)*[x(2);1]);
    p_e_y(i,:) = ([cos(x(1)) sin(x(1)) 1]*ky(:,:,i)*[x(2);1]);
    p_e_z(i,:) = ([cos(x(1)) sin(x(1)) 1]*kz(:,:,i)*[x(2);1]);
end

for i = 1:n
    for j = 1:n-i
        fx = p_e_x(i,:) - p_e_x(i+j,:);
        fy = p_e_y(i,:) - p_e_y(i+j,:);
        fz = p_e_z(i,:) - p_e_z(i+j,:);
    end
end
f = [fx;fy;fz];
end
