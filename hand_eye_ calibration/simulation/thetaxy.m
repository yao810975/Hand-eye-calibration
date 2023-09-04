%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Solve thetax, thetay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% f        :  - Circular plane normal vector 1*3 ( pay attention to its direction)
                    

% OUTPUT
% thetax   :  - The value of the radian that rotates around x
% thetay   :  - The value of the radian that rotates around y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [thetax,thetay] = thetaxy(f)


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

