%CRB克拉美罗界（一维DOA估计）---阵元数M

clc
clear

snr=10;%信噪比(dB)
source_number=1;%信号源数目
sensor_number=2:12;%阵元数目
N_x=1024;
snapshot_number=N_x;%快拍数
w=pi/4;%角频率
l=2*pi*3e8/w;%波长lamda
d=0.5*l;%阵元间距
sigma = 0.4;%噪声方差

source_doa=[50]*pi/180;%信号源的波达方向

crb = zeros(1,length(sensor_number));
for kk = 1:length(sensor_number)
    A = zeros(source_number,sensor_number(kk));%构建D行M列矩阵
    for k = 1:source_number
        A(k,:) = exp(-1i*2*pi*d*sin(source_doa(k))/l*[0:sensor_number(kk)-1]);
    end
    A = A.';%作转置运算
    
    %A对theta的求导为D
    D = zeros(sensor_number(kk),source_number);
    for k1 = 1:source_number
        for k2 = 0:sensor_number(kk)-1
            D(k2+1,k1) = (-1i*2*pi*d*cos(source_doa(k1))/l*k2*pi/180)*exp(-1i*2*pi*d*sin(source_doa(k1))/l*k2);
        end
    end
    S = sqrt(sigma*10.^(snr/10))*exp(1i*w*[0:N_x-1]);
    
    P = S*S'/N_x;%入射信号的协方差矩阵
    
    R = A*P*A' + sigma*eye(sensor_number(kk));
    
    crb(:,kk) = sqrt(sigma/(2*N_x)*inv((real(D'*(eye(sensor_number(kk))-A*inv(A'*A)*A')*D).*(P*A'*inv(R)*A*P).')));
    
end

%画图
plot(sensor_number,crb(1,:), '--rd',  'LineWidth', 1.2, 'MarkerSize', 8)
xlabel('阵元数M')
ylabel('克拉美罗界CRB')
title('克拉美罗界CRB与阵元数关系曲线')
grid on








