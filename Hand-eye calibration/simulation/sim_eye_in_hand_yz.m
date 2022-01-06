clear all

%% 参数设置
% 设置待求手眼关系Tc_e
Tc_e = Td([568,170,840])*R(1,0.5,1.5);
% 设置未知量标志点位置p_b
p_b = [10, 6.23 ,9.12];

% UR5机器人参数 
DH = [  89.2     0        pi/2
        0       -425        0
        0       -392        0
        109.3    0        pi/2
        9.475    0       -pi/2
        8.25     0           0];
%% 设定误差利用求解函数求解
% e_c_1为拟合轴线矢量数据所添加误差，e_c_2为求解数据所添加误差
for i = 0:10
    if i==0
        e_c_1 =0;
    else
        e_c_1 = unifrnd(-0.3,0.3,50,3);
    end
    e_c_2 = unifrnd(0,0.3*i,100,3);
    [e1(i+1),e2(i+1),e3(i+1),e4(i+1)] = sim_test(DH,Tc_e,p_b,e_c_1,e_c_2);
end
%% 画图
x  = [0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3.0 ];
h1=plot(x,e1,'o',...
             'LineStyle','-',...
             'Color','k',...
             'LineWidth',1.5,...
             'MarkerSize',4,...
             'MarkerFaceColor','k');
hold on;
h2=plot(x,e2,'+',...
             'LineStyle','-',...
             'Color','k',...
             'LineWidth',1.5,...
             'MarkerSize',10,...
             'MarkerFaceColor','k');
         
 hold on;
 h3=plot(x,e3,'o',...
             'LineStyle','--',...
             'Color','k',...
             'LineWidth',1,...
             'MarkerSize',4,...
             'MarkerFaceColor','k');
hold on;
h4=plot(x,e4,'+',...
             'LineStyle','--',...
             'Color','k',...
             'LineWidth',1.5,...
             'MarkerSize',10,...
             'MarkerFaceColor','k');

legend([h2 h4 h1 h3],'不进行标定件位置求解只采用y方向数据','不进行标定件位置求解只采用z方向数据','求解标定件位置只采用y方向数据','求解标定件位置只采用z方向数据',...
                                                                          'Box','off' );
set(gca,'XTick',[0 ,0.3,0.6,0.9,1.2, 1.5,1.8, 2.1, 2.4, 2.7,3.0],...
                                                                'Box','off');

xlabel('误差范围e/mm')
ylabel('标定误差均值{e_f}/mm')
%% 求解函数 
function [e1_p,e2_p,e3_p,e4_p] = sim_test(DH,Tc_e,p_b,e_c_1,e_c_2)
 %% 模拟获取含误差测量数据

% 拟合圆数据 
q = [19.4,-129,-50,-85,76,0.6]; 

for i = 1:50
    q(6) = q(6)+0.2;
    Ti = fkine(DH,q);
    P(i,:) = Tc_e*inv(Ti)*[p_b 1]';
end
P(:,4) = [];
P = P + e_c_1;
% % 求解数据(不可与拟合圆数据全部共面)

qc = unifrnd(-180,180,100,6);
for i = 1:length(qc)
    Ti = fkine(DH,qc(i,:));
    T_2(:,:,i) = Ti;
    T_1(:,:,i) = inv(Ti);
    Pc(i,:) = Tc_e*inv(Ti)*[p_b,1]';
end
Pc(:,4) = [];
Pc = Pc+e_c_2;
%% 开始求解
% Tc_b=Dc_b(Tdx,Tdy,Tdz)*Rc_b(thetax,thetay,thetaz);
% 拟合圆
[f1,c1,~] = yuan(P);
% 求解旋转参数thetax，thetay；
[thetax,thetay] = thetaxy(-f1');

%%  求解方式一只采用y方向数据
% 获取系数矩阵
k = [];
for i = 1:length(Pc)
    [ki,mi,wi,kyi] = coefficient_1(T_1(:,:,i),thetax,thetay,Pc(i,:),f1,c1);
%  只要y方向数据
    K_y(:,:,i) = mi;
    k_y(i) = kyi(2);
%  只要z方向数据
    K_z(:,:,i) = wi;
    k_z(i) = kyi(3);  
    
end

%迭代求解
% x(1) = Tdy; x(2)=p_e_x;   x(3)=p_e_z;  x(4)=p_e_y; x(5) = thetaz ;
options = optimoptions('fsolve','Display','iter');
options = optimoptions(options,'StepTolerance',1e-16);

x1 = fsolve(@(x) myfun_1(x,K_y,k_y,size(Pc)),[0;0;0;0;0],options);
p_b_1 = [x1(2),x1(4),x1(3)];
% 根据直线方程求解Tdx，Tdz
Tdy_1 = x1(1);
Tdx_1 = f1(1)*(Tdy_1-c1(2))/f1(2)+c1(1);
Tdz_1 = f1(3)*(Tdy_1-c1(2))/f1(2)+c1(3);
Tc_e_1 = Td([Tdx_1,Tdy_1,Tdz_1])*R(thetax,thetay,x1(5));
% % 求解完毕
%%  求解方式一只采用z方向数据
x3 = fsolve(@(x) myfun_1(x,K_z,k_z,size(Pc)),[0;0;0;0;0],options);
p_b_3 = [x3(2),x3(4),x3(3)];
% 根据直线方程求解Tdx，Tdz
Tdy_3 = x3(1);
Tdx_3 = f1(1)*(Tdy_3-c1(2))/f1(2)+c1(1);
Tdz_3 = f1(3)*(Tdy_3-c1(2))/f1(2)+c1(3);
Tc_e_3 = Td([Tdx_3,Tdy_3,Tdz_3])*R(thetax,thetay,x3(5));
% % 方法一只采用x方向数据
%%  求解方式二只采用y方向数据
% 
for i = 1:length(Pc)
    [kxi,kyi,kzi] = coefficient_2(T_2(:,:,i),thetax,thetay,Pc(i,:),f1,c1);
    kx(:,:,i) = kxi;
    ky(:,:,i) = kyi;
    kz(:,:,i) = kzi;
end


x2 = fsolve(@(x) myfun_2(x,ky,ky,ky,length(Pc)),[0;0],options);

Tdy_2 = x2(2);
Tdx_2 = f1(1)*(Tdy_2-c1(2))/f1(2)+c1(1);
Tdz_2 = f1(3)*(Tdy_2-c1(2))/f1(2)+c1(3);
Tc_e_2 = Td([Tdx_2,Tdy_2,Tdz_2])*R(thetax,thetay,x2(1));
%%  求解方式二只采用z方向数据
x4 = fsolve(@(x) myfun_2(x,kz,kz,kz,length(Pc)),[0;0],options);

Tdy_4 = x4(2);
Tdx_4 = f1(1)*(Tdy_4-c1(2))/f1(2)+c1(1);
Tdz_4 = f1(3)*(Tdy_4-c1(2))/f1(2)+c1(3);
Tc_e_4 = Td([Tdx_4,Tdy_4,Tdz_4])*R(thetax,thetay,x4(1));

%% 精度验证

%随机生产100个测量点验证精度
q = unifrnd(-180,180,100,6);
for i =1:100 
    Ti = fkine(DH,q(i,:)*pi/180);
    pc_test =  Tc_e*Ti*[p_b,1]';
%   获取真实Pe值
    pe = inv(Tc_e)*pc_test;
%   求解方式一标定误差
    pe_1 = inv(Tc_e_1)*pc_test;
    e1(i) = norm(pe_1 - pe);
%   求解方式一只用x方向  
    pe_3 = inv(Tc_e_3)*pc_test;
    e3(i) = norm(pe_3 - pe);
%   求解方式二标定误差    
    pe_2 = inv(Tc_e_2)*pc_test;
    e2(i)= norm(pe_2 - pe);
%   求解方式二只用x方向  
    pe_4 = inv(Tc_e_4)*pc_test;
    e4(i) = norm(pe_4 - pe);    
end
% 单位毫米
e1_p = mean(e1);
e2_p = mean(e2);
e3_p = mean(e3);
e4_p = mean(e4);
max(e1)
max(e2)
end


