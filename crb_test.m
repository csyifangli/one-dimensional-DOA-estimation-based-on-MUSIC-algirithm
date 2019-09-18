%CRB克拉美罗界（一维DOA估计）

clc
clear

snr=-10:2:10;%信噪比(dB)
source_number=1;%信号源数目
sensor_number=8;%阵元数目
N_x=1024;
snapshot_number=N_x;%快拍数
w=pi/4;%角频率
l=2*pi*3e8/w;%波长lamda
d=0.5*l;%阵元间距
sigma = 0.4;%噪声方差

source_doa=[50]*pi/180;%信号源的波达方向
A = zeros(source_number,sensor_number);%构建D行M列矩阵
for k = 1:source_number
      A(k,:) = exp(-1i*2*pi*d*sin(source_doa(k))/l*[0:sensor_number-1]);
end
A = A.';%作转置运算

%A对theta的求导为D
D = zeros(sensor_number,source_number);
for k1 = 1:source_number
    for k2 = 0:sensor_number-1
        D(k2+1,k1) = (-1i*2*pi*d*cos(source_doa(k1))/l*k2*pi/180)*exp(-1i*2*pi*d*sin(source_doa(k1))/l*k2);
    end
end

crb = zeros(1,length(snr));
for kk = 1:length(snr)
    S = sqrt(sigma*10.^(snr(kk)/10))*exp(1i*w*[0:N_x-1]);
       
    P = S*S'/N_x;%入射信号的协方差矩阵
    
    R = A*P*A' + sigma*eye(sensor_number);
    
    crb(:,kk) = sqrt(sigma/(2*N_x)*inv((real(D'*(eye(sensor_number)-A*inv(A'*A)*A')*D).*(P*A'*inv(R)*A*P).')));
    
end

%画图
plot(snr,crb(1,:), '--rd',  'LineWidth', 1.2, 'MarkerSize', 8)
% semilogy(snr,crb(1,:))
xlabel('信噪比SNR/dB')
ylabel('克拉美罗界CRB')
title('克拉美罗界CRB与信噪比关系曲线')
grid on








