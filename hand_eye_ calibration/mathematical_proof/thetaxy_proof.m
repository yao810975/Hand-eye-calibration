clear
clc
syms f1 f2 f3 real
%f_C is the Z-axis direction vector of Base in the C coordinate system
f_C = [f1;f2;f3];
%% Solve for thetax£¬thetay by f
thetax = -atan(f2/f3);
R_thetax = [    1           0               0   
          0    cos(thetax)  -sin(thetax)  
          0    sin(thetax)   cos(thetax)  ];

f_1 = inv(R_thetax)*f_C;     
thetay = atan(f_1(1)/f_1(3));

%% verify f_B = (0,0,1) = (R(thetax)*R(thetay))'*f_C
R_thetay = [cos(thetay)     0      sin(thetay)  
          0           1               0   
      -sin(thetay)    0      cos(thetay)  ];
  
simplify((R_thetax*R_thetay)'*f_C)