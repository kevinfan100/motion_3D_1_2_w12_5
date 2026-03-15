%The program is for path control

clear; clc; close all;

%parameters for simulation

Ts = 1/1600;
%lamdaC = 0.7;  % 0,4
%gammaN = 0.0425;
 

%ax_normal = 0.0147*10^6; % ( m/(N sec))

model = 'motion_3D_1_2_w12_5';

% 1. ���} Simulink �ҫ�
		load_system(model);      % ���J�ҫ��]���@�w�|���}�����^
		open_system(model);      % ���}�ҫ�����

% 2. �]�w����ɶ�
		set_param(model, 'StopTime', '3.5');
		
%---- define the ParametersBus begin----
clear elms;
elms(1) = Simulink.BusElement; elms(1).Name='Ts'; elms(1).DataType='double'; elms(1).Dimensions=1;
elms(2) = Simulink.BusElement; elms(2).Name='lamdaC'; elms(2).DataType='double'; elms(2).Dimensions=1;

elms(3) = Simulink.BusElement; elms(3).Name='theta'; elms(3).DataType='double'; elms(3).Dimensions=1;
elms(4) = Simulink.BusElement; elms(4).Name='phi'; elms(4).DataType='double'; elms(4).Dimensions=1;
elms(5) = Simulink.BusElement; elms(5).Name='pz'; elms(5).DataType='double'; elms(5).Dimensions=1;

elms(6) = Simulink.BusElement; elms(6).Name='R'; elms(6).DataType='double'; elms(6).Dimensions=1;

elms(7) = Simulink.BusElement; elms(7).Name='kb'; elms(7).DataType='double'; elms(7).Dimensions=1;
elms(8) = Simulink.BusElement; elms(8).Name='T'; elms(8).DataType='double'; elms(8).Dimensions=1;
elms(9) = Simulink.BusElement; elms(9).Name='ax_normal'; elms(9).DataType='double'; elms(9).Dimensions=1;
elms(10) = Simulink.BusElement; elms(10).Name='az_normal'; elms(10).DataType='double'; elms(10).Dimensions=1;
elms(11) = Simulink.BusElement; elms(11).Name='gammaN'; elms(11).DataType='double'; elms(11).Dimensions=1;

elms(12) = Simulink.BusElement; elms(12).Name='Avar'; elms(12).DataType='double'; elms(12).Dimensions=1;
elms(13) = Simulink.BusElement; elms(13).Name='Avar2'; elms(13).DataType='double'; elms(13).Dimensions=1;
elms(14) = Simulink.BusElement; elms(14).Name='Avar22'; elms(14).DataType='double'; elms(14).Dimensions=1;
elms(15) = Simulink.BusElement; elms(15).Name='Avar3'; elms(15).DataType='double'; elms(15).Dimensions=1;
elms(16) = Simulink.BusElement; elms(16).Name='Am_scaling'; elms(16).DataType='double'; elms(16).Dimensions=1;

elms(17) = Simulink.BusElement; elms(17).Name='beta'; elms(17).DataType='double'; elms(17).Dimensions=1;

elms(18) = Simulink.BusElement; elms(18).Name='lamdaF'; elms(18).DataType='double'; elms(18).Dimensions=1;
elms(19) = Simulink.BusElement; elms(19).Name='Pfz_11'; elms(19).DataType='double'; elms(19).Dimensions=1;
elms(20) = Simulink.BusElement; elms(20).Name='Pfz_22'; elms(20).DataType='double'; elms(20).Dimensions=1;
elms(21) = Simulink.BusElement; elms(21).Name='Pfz_33'; elms(21).DataType='double'; elms(21).Dimensions=1;
elms(22) = Simulink.BusElement; elms(22).Name='Pfz_44'; elms(22).DataType='double'; elms(22).Dimensions=1;
elms(23) = Simulink.BusElement; elms(23).Name='Pfz_55'; elms(23).DataType='double'; elms(23).Dimensions=1;
elms(24) = Simulink.BusElement; elms(24).Name='Pfz_66'; elms(24).DataType='double'; elms(24).Dimensions=1;
elms(25) = Simulink.BusElement; elms(25).Name='Pfz_77'; elms(25).DataType='double'; elms(25).Dimensions=1;

elms(26) = Simulink.BusElement; elms(26).Name='rz11_scaling'; elms(26).DataType='double'; elms(26).Dimensions=1;
elms(27) = Simulink.BusElement; elms(27).Name='rz22_scaling'; elms(27).DataType='double'; elms(27).Dimensions=1;

elms(28) = Simulink.BusElement; elms(28).Name='qz11_scaling'; elms(28).DataType='double'; elms(28).Dimensions=1;
elms(29) = Simulink.BusElement; elms(29).Name='qz22_scaling'; elms(29).DataType='double'; elms(29).Dimensions=1;
elms(30) = Simulink.BusElement; elms(30).Name='qz33_scaling'; elms(30).DataType='double'; elms(30).Dimensions=1;
elms(31) = Simulink.BusElement; elms(31).Name='qz44_scaling'; elms(31).DataType='double'; elms(31).Dimensions=1;
elms(32) = Simulink.BusElement; elms(32).Name='qz55_scaling'; elms(32).DataType='double'; elms(32).Dimensions=1;
elms(33) = Simulink.BusElement; elms(33).Name='qz66_scaling'; elms(33).DataType='double'; elms(33).Dimensions=1;
elms(34) = Simulink.BusElement; elms(34).Name='qz77_scaling'; elms(34).DataType='double'; elms(34).Dimensions=1;

ParametersBus = Simulink.Bus;
ParametersBus.Elements = elms;

assignin('base','ParametersBus',ParametersBus);

%---- define the ParametersBus end----		
		
		

% 3. �������
		out = sim(model);
		
    x_d_data = x_d;
    y_d_data = y_d;
    z_d_data = z_d;
		p_x_data = p(:,1);
		p_y_data = p(:,2);
		p_z_data = p(:,3);
		fd_x_data = squeeze(fd(1,1,:));
		fd_y_data = squeeze(fd(2,1,:));
		fd_z_data = squeeze(fd(3,1,:));

    % �ھڨ��˲v�ͦ��ɶ��b�]�P�ƾڦP�B�^
    N = size(x_d_data, 1);
    t = (0:N-1)' * Ts;

%    x_d_data = out.x_d;
%    y_d_data = out.y_d;
%    z_d_data = out.z_d;
%    fd_data = out.fd;
%    p_data = out.p;


% 4. ��X
% -----------  x ----------
clc;
figure;
subplot(311)
plot(t, x_d_data, 'o', t, p_x_data, 'r'); grid;

title('path in x '); 
xlabel('Time (s)');
ylabel('p_x (10^-6 m) ');

subplot(312)
plot(t, (x_d_data-p_x_data), 'g'); grid;

title('motion error in x '); 
xlabel('Time (s)');
ylabel('error_x (10^-6 m) ');

subplot(313)
plot(t, fd_x_data, 'r'); ; grid;
axis([0, max(t), -10, 10]);   

title('control effort in x '); 
xlabel('Time (s)');
ylabel('fd_x (pN) ');

% -----------  y ----------
clc;
figure;
subplot(311)
plot(t, y_d_data, 'o', t, p_y_data, 'r'); grid;

title('path in y '); 
xlabel('Time (s)');
ylabel('p_y (10^-6 m) ');

subplot(312)
plot(t, (y_d_data-p_y_data), 'g'); grid;

title('motion error in y '); 
xlabel('Time (s)');
ylabel('error_y (10^-6 m) ');

subplot(313)
plot(t, fd_y_data, 'r'); ; grid;
axis([0, max(t), -10, 10]);   

title('control effort in y '); 
xlabel('Time (s)');
ylabel('fd_y (pN) ');

% -----------  z ----------
clc;
figure;
subplot(311)
plot(t, z_d_data, 'o', t, p_z_data, 'r'); grid;

title('path in z '); 
xlabel('Time (s)');
ylabel('p_z (10^-6 m) ');

subplot(312)
plot( t, dz_k2, 'r'); grid;

title('motion error in z '); 
xlabel('Time (s)');
ylabel('error_z (10^-6 m) ');

subplot(313)
plot(t, fd_z_data, 'r'); ; grid;
axis([0, max(t), -10, 10]);   

title('control effort in z '); 
xlabel('Time (s)');
ylabel('fd_z (pN) ');

% -----------  z motion gain ----------

clc;
figure;
plot(t, azm_k, 'o', t, az_hat_k, 'r', t, mgain_z, 'g'); grid;
%axis([0, max(t), -1.1*max(abs(azm_k)*10^-6), 1.1*max(abs(azm_k)*10^-6)]);   

title('path in z '); 
xlabel('Time (s)');
ylabel('az (um/(pN s)) ');

clc;
figure;

%plot( t, 1.4*az_hat_k, 'r', t, mgain_z, 'g'); grid;
plot( t, 1.0*az_hat_k, 'r', t, mgain_z, 'g'); grid;


title('path in z '); 
xlabel('Time (s)');
ylabel('az (um/(pN s)) ');


% -----------  estimator feedback matrix L ----------
 
clc;
figure;
subplot(231)
plot(t, squeeze(Lz(1,1,:)), 'o', t, squeeze(Lz(2,1,:)), 'r', t, squeeze(Lz(3,1,:)), 'g'); grid;
%plot(t, Lz(1,1)), 'o', t, Lz(2,1), 'r', t, Lz(3,1), 'g'); grid;
%axis([0, max(t), -0.1, 1.1]);   

title('path in z '); 
xlabel('Time (s)');
ylabel('Lz1 '); 
 
subplot(234)
plot(t, squeeze(Lz(1,2,:)), 'o', t, squeeze(Lz(2,2,:)), 'r', t, squeeze(Lz(3,2,:)), 'g'); grid;
%axis([0, max(t), -0.1, 1.1]);   

title('path in z '); 
xlabel('Time (s)');
ylabel('Lz2 '); 

subplot(232)
plot(t, squeeze(Lz(4,1,:)), 'o', t, squeeze(Lz(5,1,:)), 'r'); grid;
%axis([0, max(t), -0.1, 1.1]);   

title('path in z '); 
xlabel('Time (s)');
ylabel('Lzd1 '); 
 
subplot(235)
plot(t, squeeze(Lz(4,2,:)), 'o', t, squeeze(Lz(5,2,:)), 'r'); grid;
%axis([0, max(t), -0.1, 1.1]);   

title('path in z '); 
xlabel('Time (s)');
ylabel('Lzd2 ');  

subplot(233)
plot(t, squeeze(Lz(6,1,:)), 'o', t, squeeze(Lz(7,1,:)), 'r'); grid;
%axis([0, max(t), -0.1, 1.1]);   

title('path in z '); 
xlabel('Time (s)');
ylabel('Laz1 '); 
 
subplot(236)
plot(t, squeeze(Lz(6,2,:)), 'o', t, squeeze(Lz(7,2,:)), 'r'); grid;
%axis([0, max(t), -0.1, 1.1]);   

title('path in z '); 
xlabel('Time (s)');
ylabel('Laz2 ');  


% -----------  forecast of error (state) covariance matrix Pfz ----------
 
clc;
figure;
subplot(311)
plot(t, squeeze(Pfz_k(1,1,:)), 'o', t, squeeze(Pfz_k(2,2,:)), 'r', t, squeeze(Pfz_k(3,3,:)), 'g'); grid;
%axis([0, max(t), -0.1, 1.1]);   

title('path in z '); 
xlabel('Time (s)');
ylabel('Pfz_dx '); 

subplot(312)
plot(t, squeeze(Pfz_k(4,4,:)), 'o', t, squeeze(Pfz_k(5,5,:)), 'r'); grid;
%axis([0, max(t), -0.1, 1.1]);   

title('path in z '); 
xlabel('Time (s)');
ylabel('Pfz_xd '); 

subplot(313)
plot(t, squeeze(Pfz_k(6,6,:)), 'o', t, squeeze(Pfz_k(7,7,:)), 'r'); grid;
%axis([0, max(t), -0.1, 1.1]);   

title('path in z '); 
xlabel('Time (s)');
ylabel('Pfz_az '); 

% -----------  estimated state  ----------
 
clc;
figure;
subplot(311)
plot( t, dz1_hat_k, 'r'); grid;

title('motion error in z '); 
xlabel('Time (s)');
ylabel('dz1_hat_k (10^-6 m) ');

subplot(312)
plot( t, dz2_hat_k, 'r'); grid;

title('motion error in z '); 
xlabel('Time (s)');
ylabel('dz2_hat_k (10^-6 m) ');

subplot(313)
plot( t, dz3_hat_k, 'r'); grid;

title('motion error in z '); 
xlabel('Time (s)');
ylabel('dz3_hat_k (10^-6 m) ');

clc;
figure;
plot( t, dz1_hat_k, 'b', t, dz2_hat_k,'r', t, dz3_hat_k,'g'); grid;
title('motion error in z ');
xlabel('Time (s)');
ylabel('dz_hat_k (10^-6 m) ');

% -----------  3D trajectory plot (moved from Simulink plot3_block) ----------
figure;
plot3(p_x_data, p_y_data, p_z_data, '.-'); hold on;
plot3(x_d_data, y_d_data, z_d_data, 'r--');
grid on; axis equal;
xlabel('x (um)'); ylabel('y (um)'); zlabel('z (um)');
title('3D Trajectory');
legend('actual', 'desired');
view(3); rotate3d on;
