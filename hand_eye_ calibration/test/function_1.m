%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  equation set£ºf =  [cos(x(5)),sin(x(5)),1]*Kx*[x(2);x(3);x(4);1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% x        :  - Unknown parameter£ºx(1) = y1; x(2)=x2;  x(3)=z2;  x(4)=y2;  x(5) = thetaz ;
% K        :  - Coefficient matrix (excluding y1)
% k        :  - Coefficient matrix of y1
% n        :  - The number of equations n=3*size(P);

% OUTPUT
% f        :  - The system to be solved by 'fsolve' function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=function_1(x,K,k,n)

Kn =[]; 
for i = 1:n
    Kx = K(:,:,i);
    Kx(3,4) = k(i)*x(1) + Kx(3,4);
    Kwi = [cos(x(5)),sin(x(5)),1]*Kx;
    Kn = [Kn;Kwi];
end

 f= Kn*[x(2);x(3);x(4);1];

end     
