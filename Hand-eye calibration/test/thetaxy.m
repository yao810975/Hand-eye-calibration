function [thetax,thetay] = thetaxy(f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 代码采用先绕x后绕y的几何法求解R(thetax,thetay,thetaz)中的thetax，thetay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% b        :  - 圆平面法向量,(输入方向需为旋转轴的z轴正方向)

% OUTPUT
% thetax   :  - 绕x旋转的弧度值
% thetay   :  - 绕y旋转的弧度值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if f(3)<0
    thetax = -atan(f(2)/f(3))+pi;

else
    thetax = -atan(f(2)/f(3));

end

Rx = [    1           0               0   
          0    cos(thetax)  -sin(thetax)  
          0    sin(thetax)   cos(thetax)  ];
fx = inv(Rx)*f';
if fx(3)<0
    
    thetay = atan(fx(1)/fx(3))+pi;
else

    thetay = atan(fx(1)/fx(3));
end

end

