function [R] = R(thetax,thetay,thetaz)
%旋转变换方程
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

