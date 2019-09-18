%阵元间距对MUSIC算法的影响
clc
clear 
format long
lamda = 150;%最高频率信号的波长

d = lamda/6;%阵元间距
d2 = lamda/2;
d3 = lamda;
%此处需注意length(theta) = length(w)
theta = [20,30,60]/180*pi;%信号入射角度
w = [pi/6,pi/4,pi/3];%角频率
w = w';
snapshots = 200;%快拍数
D = length(w);%信号源数目
M = 12;%天线阵元数目
SNR = 20;%信噪比为20dB
A = zeros(D,M);%构建D行M列矩阵
for k = 1:D
      A(k,:) = exp(-1i*2*pi*d*sin(theta(k))/lamda*[0:M-1]);
      A2(k,:) = exp(-1i*2*pi*d2*sin(theta(k))/lamda*[0:M-1]);
      A3(k,:) = exp(-1i*2*pi*d3*sin(theta(k))/lamda*[0:M-1]);
end

A = A.';%获得方向矩阵
A2 = A2.';
A3 = A3.';
S = 4*exp(1i*(w*[1:snapshots]));%获得仿真信号(根据窄带信号的定义式)
X = awgn(A*S,SNR,'measured');%接收信号
X2 = awgn(A2*S,SNR,'measured');
X3 = awgn(A3*S,SNR,'measured');
Rx = X*X'/snapshots;%接收信号的协方差矩阵
Rx2 = X2*X2'/snapshots;
Rx3 = X3*X3'/snapshots;
[Ve,Va] = eig(Rx);%获取协方差矩阵的特征值Va以及特征向量Ve
En = Ve(:,1:M-D);%取出M-D个对应特征值为噪声方差的特征向量

[Ve2,Va2] = eig(Rx2);%获取协方差矩阵的特征值Va以及特征向量Ve
En2 = Ve2(:,1:M-D);%取出M-D个对应特征值为噪声方差的特征向量

[Ve3,Va3] = eig(Rx3);%获取协方差矩阵的特征值Va以及特征向量Ve
En3 = Ve3(:,1:M-D);%取出M-D个对应特征值为噪声方差的特征向量


theta1 = -90:0.5:90;


%使用两个for循环进行嵌套，将各个角度依次带入进行谱峰搜索
for a = 1:length(theta1)
    AA = zeros(1,M);
    AA2 = zeros(1,M);
    AA3 = zeros(1,M);
    for b = 0:M-1
        AA(:,b+1) = exp(-1*1i*2*pi*d*sin(theta1(a)/180*pi)/lamda*b); %注：在matlab中三角函数计算使用地是弧度值
         AA2(:,b+1) = exp(-1*1i*2*pi*d2*sin(theta1(a)/180*pi)/lamda*b);
          AA3(:,b+1) = exp(-1*1i*2*pi*d3*sin(theta1(a)/180*pi)/lamda*b);
    end
    AA = AA.';
    AA2 = AA2.';
    AA3 = AA3.';
    
    P = AA'*En*En'*AA;
    P2 = AA2'*En2*En2'*AA2;
    P3 = AA3'*En3*En3'*AA3;
    Pmusic(a) = abs(1/P);%谱函数
    Pmusic2(a) = abs(1/P2);
    Pmusic3(a) = abs(1/P3);
end

Pmusic = 10*log10(Pmusic/max(Pmusic));%将获得的谱函数归一化后进行计算dB操作
Pmusic2 = 10*log10(Pmusic2/max(Pmusic2));
Pmusic3 = 10*log10(Pmusic3/max(Pmusic3));
plot(theta1,Pmusic,'-r')
hold on
plot(theta1,Pmusic2,'--g')
hold on
plot(theta1,Pmusic3,'-.b')
hold off

xlabel('角度 \theta/degree')
ylabel('谱函数 P(\theta)/dB')
title('基于MUSIC算法的DOA估计')
legend('阵元间距为150/6','阵元间距为150/2','阵元间距为150')
text(-90,-20,'波长为150')
grid on






