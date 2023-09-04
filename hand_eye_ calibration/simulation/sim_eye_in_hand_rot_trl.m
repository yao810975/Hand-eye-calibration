clear 
clc
%% Parameter settings
%  set the Tc_e of the desired hand-eye transformation and the position of the marker p_b

wc_e = [1,0.5,1.5];
tc_e = [10 -19 107];
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

    for k = 1:1000
   %% Simulation acquires measurement data with errors
           
        e_c = normrnd(0,0.2*j,100,3);
        e_cir = normrnd(0,0.2*j,100,3);
        
        % Fit circle data (set to joint 6 from largest to small rotation
        % 300бу)
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
        %% Bring in the data to solve
         % solution 1 and solution 2         
        [Tc_e_1,Tc_e_2,f,c]=eye_in_hand(DH,Pcir,Pcal,qc);
         % conventional method         
        [Tc_e_3] = Kronecker_Product(DH,Pcal,qc);

        % correct rotation part of Tc_e_3 
        Rc_e_3 = Tc_e_3(1:3,1:3);
        [wc_e_3(1),wc_e_3(2),wc_e_3(3)] = dcm2angle(Rc_e_3','XYZ');
        [thetax,thetay] = thetaxy(f');
        wc_e_4 = [thetax,thetay,wc_e_3(3)];  
        Rc_e_4 = angle2dcm(wc_e_4(1),wc_e_4(2),wc_e_4(3),'XYZ')';
        Tc_e_4 = [Rc_e_4; [0 0 0]];
        Tc_e_4 = [Tc_e_4 Tc_e_3(:,4)];
        
        %% Accuracy verification
        [e1_rot(k),e1_t(k)] = error_test(wc_e,tc_e,Tc_e_1);  
        [e2_rot(k),e2_t(k)] = error_test(wc_e,tc_e,Tc_e_2);
        [e3_rot(k),e3_t(k)] = error_test(wc_e,tc_e,Tc_e_3);
        [e4_rot(k),e4_t(k)] = error_test(wc_e,tc_e,Tc_e_4);
    end
    
    e1_rot_m(j+1) = mean(e1_rot);
    e1_t_m(j+1) = mean(e1_t);
    e2_rot_m(j+1) = mean(e2_rot);
    e2_t_m(j+1) = mean(e2_t);
    e3_rot_m(j+1) = mean(e3_rot);
    e3_t_m(j+1) = mean(e3_t);
    e4_rot_m(j+1) = mean(e4_rot);
    e4_t_m(j+1) = mean(e4_t);
end
%% plot
figure(1)
x  = [0,1,2,3,4,5,6,7,8,9,10];
h1=plot(x,e1_rot_m,'o',...
             'LineStyle','-',...
             'Color','k',...
             'LineWidth',1.5,...
             'MarkerSize',4,...
             'MarkerFaceColor','k');
hold on;
h2=plot(x,e2_rot_m,'+',...
             'LineStyle','-',...
             'Color','k',...
             'LineWidth',1.5,...
             'MarkerSize',10,...
             'MarkerFaceColor','k');
h3=plot(x,e3_rot_m,'o',...
             'LineStyle','--',...
             'Color','b',...
             'LineWidth',1.5,...
             'MarkerSize',4,...
             'MarkerFaceColor','b');                 
h4=plot(x,e4_rot_m,'+',...
             'LineStyle','-',...
             'Color','b',...
             'LineWidth',1.5,...
             'MarkerSize',10,...
             'MarkerFaceColor','b'); 

legend([h3 h4 h2 h1],'Kronecker Product','corrected by \thetax\thetay','solution 1','solution 2',...
                                                     'Box','off' );
set(gca,'XTick',[0 ,1,2,3,4, 5,6, 7, 8, 9, 10]);

xlabel('noise level')
ylabel('rotation error(rad)')

figure(2)
x  = [0,1,2,3,4,5,6,7,8,9,10];
h1=plot(x,e1_t_m,'o',...
             'LineStyle','-',...
             'Color','k',...
             'LineWidth',1.5,...
             'MarkerSize',4,...
             'MarkerFaceColor','k');
hold on;
h2=plot(x,e2_t_m,'+',...
             'LineStyle','-',...
             'Color','k',...
             'LineWidth',1.5,...
             'MarkerSize',10,...
             'MarkerFaceColor','k');
h3=plot(x,e3_t_m,'o',...
             'LineStyle','--',...
             'Color','b',...
             'LineWidth',1.5,...
             'MarkerSize',4,...
             'MarkerFaceColor','b');                 
h4=plot(x,e4_t_m,'+',...
             'LineStyle','-',...
             'Color','b',...
             'LineWidth',1.5,...
             'MarkerSize',10,...
             'MarkerFaceColor','b'); 

legend([h3 h4 h2 h1],'Kronecker Product','corrected by \thetax\thetay','solution 1','solution2',...
                                                     'Box','off' );
set(gca,'XTick',[0 ,1,2,3,4, 5,6, 7, 8, 9, 10]);

xlabel('noise level')
ylabel('translation error(mm)')

% Accuracy verification
function [e_rot,e_t] = error_test(wc_e,tc_e,Tc_e_x)
%    rotation error 
    Rc_e_x = Tc_e_x(1:3,1:3);
    [wc_e_x(1),wc_e_x(2),wc_e_x(3)] = dcm2angle(Rc_e_x','XYZ');
%    translation error
    tc_e_x = Tc_e_x(1:3,4);
    e_rot = norm(wc_e_x - wc_e);
    e_t = norm(tc_e' - tc_e_x);
 
end


