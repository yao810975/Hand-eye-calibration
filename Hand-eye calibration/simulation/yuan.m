function [A,C,r] = yuan(Data)

M = Data;
[num,dim]=size(M);
 
L1=ones(num,1);
A=inv(M'*M)*M'*L1       % 求解平面法向量

B=zeros((num-1)*num/2,3);
 
count=0;
for i=1:num-1
    for j=i+1:num   
        count=count+1;
        B(count,:)=M(j,:)-M(i,:);
    end    
end
 
L2=zeros((num-1)*num/2,1);
count=0;
for i=1:num-1
    for j=i+1:num
        count=count+1;
        L2(count)=(M(j,1)^2+M(j,2)^2+M(j,3)^2-M(i,1)^2-M(i,2)^2-M(i,3)^2)/2;
    end
end
 
D=zeros(4,4);
D(1:3,1:3)=(B'*B);
D(4,1:3)=A';
D(1:3,4)=A;
 
L3=[B'*L2;1];
 
C=inv(D')*(L3);   % 求解空间圆圆心坐标
 
C=C(1:3);
 
radius=0;
for i=1:num
    tmp=M(i,:)-C';
    radius=radius+sqrt(tmp(1)^2+tmp(2)^2+tmp(3)^2);
end
r=radius/num;            %  空间圆拟合半径
 
 
 


end

