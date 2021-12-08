function [T06] = fkine(DH,q)
% 正运动学
d = DH(:,1);
a = DH(:,2);
alpha = DH(:,3);
T01 =mT(q(1),  d(1) ,    a(1),     alpha(1));%[theta d a alpha ]  
T12 =mT(q(2),  d(2),     a(2),     alpha(2));
T23 =mT(q(3),  d(3) ,    a(3),     alpha(3));
T34 =mT(q(4),  d(4) ,    a(4),     alpha(4));
T45 =mT(q(5),  d(5) ,    a(5),     alpha(5));
T56 =mT(q(6),  d(6) ,    a(6),     alpha(6));
T06=T01*T12*T23*T34*T45*T56;

function [T] = mT(theta,d,a,aplha)
% 机器人正运动学公式
Rx = [1         0               0             0
    0       cos(aplha)       -sin(aplha)      0
    0       sin(aplha)        cos(aplha)      0
    0         0                 0             1];

  
Rz = [cos(theta)    -sin(theta)     0   0
      sin(theta)    cos(theta)      0   0
          0         0               1   0
          0         0               0   1];

Dx  = [1        0         0         a
      0         1         0        0
      0         0         1        0
      0         0         0        1];

Dz  = [1        0         0         0
      0         1         0        0
      0         0         1        d
      0         0         0        1];
% MDH  
% T = Rx*Dx*Rz*Dz;
% SDH
T = Dz*Rz*Dx*Rx;
end


end

