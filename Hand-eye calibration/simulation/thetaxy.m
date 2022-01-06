function [thetax,thetay] = thetaxy(f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % �����������x����y�ļ��η����R(thetax,thetay,thetaz)�е�thetax��thetay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% b        :  - Բƽ�淨���� 1*3(����Ϊ����z��������,��⵼�������󣬿��Ǵ˴���������

% OUTPUT
% thetax   :  - ��x��ת�Ļ���ֵ
% thetay   :  - ��y��ת�Ļ���ֵ
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

