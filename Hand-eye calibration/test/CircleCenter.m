%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��Բ������ȷ��Բ�ĺͰ뾶
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% p1   :  - ��һ��������, ������ 1x3
% p2   :  - �ڶ���������, ������ 1x3
% p3   :  - ������������, ������ 1x3
% OUTPUT
% pc   :  - Բ������, ������ 1x3
% r    :  - �뾶, ����
% f    :  - Բƽ�淨����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [pc,r,f]=CircleCenter(p1,p2,p3)

    % ������
    validateattributes(p1,{'numeric'},{'row'},1);% ������
    validateattributes(p2,{'numeric'},{'row'},2);
    validateattributes(p3,{'numeric'},{'row'},3);
    num1=length(p1);num2=length(p2);num3=length(p3);
    if (num1 == num2) && (num2 == num3)
        if num1 == 2
            p1=[p1,0];p2=[p2,0];p3=[p3,0];
        elseif num1 ~= 3
            error('��֧�ֶ�ά����ά��������');
        end
    else
        error('���������ά����һ��');
    end

    % ���߼��
    temp01=p1-p2;temp02=p3-p2;
    temp03=cross(temp01,temp02);
    temp=(temp03*temp03')/(temp01*temp01')/(temp02*temp02');
    if temp < 10^-5
        pc = 0; r=0; f=[0 0 0]';
    end
    if temp < 10^-10
        error('���㹲��, �޷�ȷ��Բ');
    end

    mat1=[p1,1;p2,1;p3,1];% size = 3x4

    m=+det(mat1(:,2:4));
    n=-det([mat1(:,1),mat1(:,3:4)]);
    p=+det([mat1(:,1:2),mat1(:,4)]);
    q=-det(mat1(:,1:3));

    mat2=[[p1*p1';p2*p2';p3*p3'],mat1;2*q,[-m,-n,-p,0]];% size = 4x5

    A=+det(mat2(:,2:5));
    B=-det([mat2(:,1),mat2(:,3:5)]);
    C=+det([mat2(:,1:2),mat2(:,4:5)]);
    D=-det([mat2(:,1:3),mat2(:,5)]);
    E=+det(mat2(:,1:4));

    pc =-[B,C,D]/2/A;
    r=sqrt(B^2+C^2+D^2-4*A*E)/2/abs(A);
    f = [m n p];
    
end
