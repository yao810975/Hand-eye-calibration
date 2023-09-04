%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Kronecker Product method


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% Pcam            :  - Measurement 3D data
% DH and qc       :  - Corresponding robot pose (from end to base)

% OUTPUT
% Tc_e   £º - Hand-eye calibration result

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Tc_e] = Kronecker_Product(DH,Pcam,qc)
    
   for i = 1:size(qc,1)
       
       Tb_e_i = fkine(DH,qc(i,:));
       Te_b_i = inv(Tb_e_i);
       Te_b(:,:,i) = Te_b_i;
   end
  
    K =[];
    b =[];
    for i = 1:size(Pcam,1)
        
        Te_b_i = Te_b(:,:,i);
        Re_b_i = Te_b_i(1:3,1:3);
        De_b_i = Te_b_i(1:3,4);
        k_i = [kron(Pcam(i,:),eye(3,3)), eye(3,3), -Re_b_i];
        K= [K;k_i];
        b_i = De_b_i;
        b = [b;b_i];
    end
    x = inv((K'*K))*K'*b;

    Re_cx = reshape(x(1:9),3,3);
    [U,S,V] = svd(Re_cx);
    Re_cx = sign(det(Re_cx))*U*V';
    te_c = x(10:12);
    
    Re_c = [Re_cx;[0 0 0]];
    Te_c = [Re_c [te_c;1]] ;
    Tc_e = inv(Te_c);
end

