
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Fit the circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% Data     :  - Fit circle data £¨n*3£©

% OUTPUT
% A        :  - Circular plane normal vector
% C        :  - Center
% r        :  - radius
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [normal,C,r] = fit_circle(Data)

M = Data;
[num,dim]=size(M);
L1=ones(num,1);
A=inv(M'*M)*M'*L1;       % Solve for plane normal vectors

B=zeros((num-1)*num/2,3);
 
count=0;
for i=1:num-1
    for j=i+1:num   
        count=count+1;
        B(count,:)=M(j,:)-M(i,:);
        a(count,:) = cross(B(i,:), B(j,:));
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
 
C=inv(D')*(L3);   % Solve for the coordinates of the center of the  circle
 
C=C(1:3);
 
radius=0;
for i=1:num
    tmp=M(i,:)-C';
    radius=radius+sqrt(tmp(1)^2+tmp(2)^2+tmp(3)^2);
end
r=radius/num;            %  Space circle fitting radius

% To remove the center of mass
W = M - repmat(mean(M),size(M,1),1);
% SVD(singular value decomposition)
[~,~,V] = svd(W,0);
normal = V(:,3);
% Look for the positive z-axis
g = median(a(:,1)); 
if sign(g)~=sign(normal(1))
    normal = -normal;
end 
 


end

