%当阵元间距大于半波长时
%MUSIC
clc
clear 
format long
lamda = 150;%最高频率信号的波长
d = 0.563*lamda;%阵元间距
%此处需注意length(theta) = length(w)
theta = [20]/180*pi;%信号入射角度
w = [pi/4];%角频率
w = w';
snapshots = 200;%快拍数
D = length(w);%信号源数目
M = 12;%天线阵元数目
SNR = 10;%信噪比为10dB
A = zeros(D,M);%构建D行M列矩阵
for k = 1:D
      A(k,:) = exp(-1i*2*pi*d*sin(theta(k))/lamda*[0:M-1]);
end

A = A.';%获得方向矩阵

S = 4*exp(1i*(w*[1:snapshots]));%获得仿真信号(根据窄带信号的定义式)
X = awgn(A*S,SNR,'measured');%接收信号
Rx = X*X'/snapshots;%接收信号的协方差矩阵
[Ve,Va] = eig(Rx);%获取协方差矩阵的特征值Va以及特征向量Ve
En = Ve(:,1:M-D);%取出M-D个对应特征值为噪声方差的特征向量

theta1 = -90:0.5:90;


%使用两个for循环进行嵌套，将各个角度依次带入进行谱峰搜索
for a = 1:length(theta1)
    AA = zeros(1,M);
    for b = 0:M-1
        AA(:,b+1) = exp(-1*1i*2*pi*d*sin(theta1(a)/180*pi)/lamda*b); %注：在matlab中三角函数计算使用地是弧度值
    end
   AA = AA.';
    P = AA'*En*En'*AA;
    Pmusic(a) = abs(1/P);%谱函数
end

Pmusic = 10*log10(Pmusic/max(Pmusic));%将获得的谱函数归一化后进行计算dB操作
plot(theta1,Pmusic,'-r')
xlabel('角度 \theta/degree')
ylabel('谱函数 P(\theta)/dB')
title('基于MUSIC算法的DOA估计')
grid on






