clear 
clc
%% Parameter settings
%  Arbitrarily set the Tc_e of the desired hand-eye transition and the
%  position of the marker point p_b

Tc_e = Td([10 -19 107])*R(1,0.5,1.5);
p_b = [568,170,440];

% DH parameters of UR5 
DH = [  89.2     0        pi/2
        0       -425        0
        0       -392        0
        109.3    0        pi/2
        90.475    0       -pi/2
        80.25     0           0];
    
    
%% Simulation
% j is the noise scale
% k is the number of calibrations

for j = 0:10
    for k =1:1000
        %% Simulation acquires measurement data with errors      

        e_c = normrnd(0,0.2*j,100,3);
        e_cir = normrnd(0,0.2*j,100,3);

        % Fit circle data (set to joint 6 from largest to small rotation
%         300бу)
        q = unifrnd(-pi,pi,1,6);

        Pcir = [];
        for i = 1:100
            q(6) = q(6)-3*pi/180;
            Ti = fkine(DH,q);
            Pcir(i,:) = Tc_e*inv(Ti)*[p_b 1]';
        end
        Pcir(:,4) = [];
        [f1,c1,~] = fit_circle(Pcir);
        Pcir = Pcir + e_cir;

        % Solve the data (100 random movements)
        Pcal = [];
        qc = unifrnd(-pi,pi,100,6);
        for i = 1:length(qc)
            Tb_e_i = fkine(DH,qc(i,:));
            Te_b_i = inv(Tb_e_i);
             % The rotation of the Te_b is converted to Euler angles and 
             %  noise simulations are added to obtain the actual value
            Re_b_i = Te_b_i(1:3,1:3);
            [omiga(1),omiga(2),omiga(3)] = dcm2angle(Re_b_i','XYZ');
            omiga = omiga + normrnd(0,0.2*j*pi/180,1,3);
            rRe_b_i = angle2dcm(omiga(1),omiga(2),omiga(3),'XYZ')';
            % Te_b Translation adds noise to the simulated actual value    
            rDe_b_i = Te_b_i(1:3,4) + normrnd(0,0.2*j,1,3)';
     
            % Combined translation and rotation simulation to obtain the 
            % actual measured robot pose
            rTe_b_i = [rRe_b_i;[0 0 0 ]];
            rTe_b_i = [rTe_b_i, [rDe_b_i;1]];       
            Pcal(i,:) = Tc_e*rTe_b_i*[p_b,1]';       
        end
        Pcal(:,4) = [];
        Pcal = Pcal+e_c;


        %% Under the same conditions, only x, y, z single direction data is used
        [Tc_e_1_x,Tc_e_1_y,Tc_e_1_z,Tc_e_2_x,Tc_e_2_y,Tc_e_2_z] = eye_in_hand_xyz(DH,Pcir,Pcal,qc);

        %% Accuracy verification
        %  Randomly generated 100 points verification accuracy (100 points)
        q_test = unifrnd(-pi,pi,100,6);
        for i =1:100 
            Ti = fkine(DH,q_test(i,:));
            Ti_all(:,:,i) = Ti;
            Pc_test(:,i) = Tc_e*inv(Ti)*[p_b 1]';
        end
        
        %  solution 1 uses only a single direction
        e1(k) = pre_test(p_b,Pc_test,Tc_e_1_x,Ti_all);
        e2(k) = pre_test(p_b,Pc_test,Tc_e_1_y,Ti_all);
        e3(k) = pre_test(p_b,Pc_test,Tc_e_1_z,Ti_all);
        %  solution 2 uses only a single direction
        e4(k) = pre_test(p_b,Pc_test,Tc_e_1_x,Ti_all);
        e5(k) = pre_test(p_b,Pc_test,Tc_e_1_y,Ti_all);
        e6(k) = pre_test(p_b,Pc_test,Tc_e_1_z,Ti_all);
    end
    
    e1_m(j+1) = mean(e1);
    e2_m(j+1) = mean(e2);
    e3_m(j+1) = mean(e3);
    e4_m(j+1) = mean(e4);
    e5_m(j+1) = mean(e5);
    e6_m(j+1) = mean(e6);

end
%% plot
x  = [0,1,2,3,4,5,6,7,8,9,10];
figure
h1=plot(x,e1_m,'o',...
             'LineStyle','-',...
             'Color','r',...
             'LineWidth',1,...
             'MarkerSize',4,...
             'MarkerFaceColor','r');
hold on;
h2=plot(x,e2_m,'o',...
             'LineStyle','-',...
             'Color','b',...
             'LineWidth',1,...
             'MarkerSize',4,...
             'MarkerFaceColor','b');
h3=plot(x,e3_m,'o',...
             'LineStyle','-',...
             'Color','k',...
             'LineWidth',1,...
             'MarkerSize',4,...
             'MarkerFaceColor','k');

h4=plot(x,e4_m,'o',...
             'LineStyle','--',...
             'Color','r',...
             'LineWidth',1.5,...
             'MarkerSize',4,...
             'MarkerFaceColor','r');
 h5=plot(x,e5_m,'o',...
             'LineStyle','--',...
             'Color','b',...
             'LineWidth',1.5,...
             'MarkerSize',4,...
             'MarkerFaceColor','b');

h6=plot(x,e6_m,'o',...
             'LineStyle','--',...
             'Color','k',...
             'LineWidth',1.5,...
             'MarkerSize',4,...
             'MarkerFaceColor','k');        


legend([h1  h2 h3 h4 h5 h6 ],'solution 1-x','solution 1-y','solution 1-z','solution 2-x','solution 2-y','solution 2-z',...
                                                          'Box','off' );
set(gca,'XTick',[0,1, 2, 3, 4, 5, 6, 7, 8,9, 10],...
                                                                'Box','off');
xlabel('noise level')
ylabel('position error(mm)')

% Accuracy verification
% position error
function [e_1] = pre_test(p_b,Pc_test,Tc_e_x,Ti_all)
    
    p_b = [p_b 1]';
    for i =1:size(Pc_test,2)
        
        Ti = Ti_all(:,:,i);
        p_b_x = Ti*inv(Tc_e_x)*Pc_test(:,i);
        
        e = norm(p_b - p_b_x);
    end
    e_1 = mean(e);
end