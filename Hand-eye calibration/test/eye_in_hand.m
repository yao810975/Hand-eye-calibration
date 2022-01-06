clear all


%% 参数设置
% 设置待求手眼关系Tc_b
Tc_e =  Td([56,17,15])*R(1,0.5,1.1);
% 设置待求标志点位置p_e
p_b = [30, 32 ,510];

% UR5机器人参数 
DH = [  89.2     0        pi/2
        0       -425        0
        0       -392        0
        109.3    0        pi/2
        9.475    0       -pi/2
        8.25     0           0];
   
    
%% 模拟获取无偏测量数据
% p1 p2 p3 设置为关节6转动用以确定拟合轴线矢量。（从大到小转动确定矢量方向）
% pc 不与p1 p2 p3共面用以最终求解
q1 = [2   -129.16 -25.67 -85.39 0 86.66];
q2 = [2  -129.16 -25.67 -85.39 0 56.66];
q3 = [2  -129.16 -25.67 -85.39 0 6.66];
qc = [0   -129.16 -25.67 -85.39 0 86.66];
[T1] = fkine(DH,q1*pi/180);      
[T2] = fkine(DH,q2*pi/180);      
[T3] = fkine(DH,q3*pi/180); 
[Tc] = fkine(DH,qc*pi/180);
% p_c = Tc_e*Te_b*p_e 
p1 = Tc_e*inv(T1)*[p_b,1]';p2 = Tc_e*inv(T2)*[p_b,1]';
p3 = Tc_e*inv(T3)*[p_b,1]';pc = Tc_e*inv(Tc)*[p_b,1]';
p1 = [p1(1) p1(2) p1(3)]; p2 = [p2(1) p2(2) p2(3)]; 
p3 = [p3(1) p3(2) p3(3)]; pc = [pc(1) pc(2) pc(3)];

P  = [p2;pc];
T(:,:,1) = T2;
T(:,:,2) = Tc;

%% 开始求解
% Tc_b=Dc_b(Tdx,Tdy,Tdz)*Rc_b(thetax,thetay,thetaz);
% 拟合圆
[c1,r1,f1] = CircleCenter(p1,p2,p3);
% 求解旋转参数thetax，thetay；   
[thetax,thetay] = thetaxy(f1);

% 求解方式一：求出标定点位置信息
% 获取K矩阵

k = [];
for i = 1:size(P,1)
    [ki,mi,wi,kyi] = coefficient_1(inv(T(:,:,i)),thetax,thetay,P(i,:),f1,c1);
    K(:,:,3*i-2) = ki;
    K(:,:,3*i-1) = mi;
    K(:,:,3*i) = wi;
    k = [k;kyi];
end

% x(1) = Tdy; x(2)=p_e_x;   x(3)=p_e_z;  x(4)=p_e_y; x(5) = thetaz ;

x = fsolve(@(x) myfun_1(x,K,k,3*size(P,1)),[0;0;0;0;0],optimoptions('fsolve','Display','iter'));

p_b_1 = [x(2),x(4),x(3)];
% 根据直线方程求解Tdx，Tdz
Tdy = x(1);
Tdx = f1(1)*(Tdy-c1(2))/f1(2)+c1(1);
Tdz = f1(3)*(Tdy-c1(2))/f1(2)+c1(3);
Tc_e_1 = Td([Tdx,Tdy,Tdz])*R(thetax,thetay,x(5));
%% 求解方式二：不计算标志点位置信息

for i = 1:size(P,1)
    [kxi,kyi,kzi] = coefficient_2(T(:,:,i),thetax,thetay,P(i,:),f1,c1);
    kx(:,:,i) = kxi;
    ky(:,:,i) = kyi;
    kz(:,:,i) = kzi;
end

%  x2(1) = thetaz ; x2(2) = Tdy; 
x2 = fsolve(@(x2) myfun_2(x2,kx,ky,kz,size(P,1)),[0;0],optimoptions('fsolve','Display','iter'));
Tdy = x2(2);
Tdx = f1(1)*(Tdy-c1(2))/f1(2)+c1(1);
Tdz = f1(3)*(Tdy-c1(2))/f1(2)+c1(3);
Tc_b_2 = Td([Tdx,Tdy,Tdz])*R(thetax,thetay,x2(1));
% 方法二求解完毕
%% 验证
p_b_1 - p_b
Tc_b_2 - Tc_e
Tc_e_1 - Tc_e








