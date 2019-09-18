%MUSIC
clc
clear 
format long
lamda = 150;%最高频率信号的波长
m1 = 5;
m2 = 6;%保证m1和m2两个数互为质数
d1 = m1*lamda/2;%阵元间距
d2 = m2*lamda/2;%阵元间距
%此处需注意length(theta) = length(w)
theta = [20,60]/180*pi;%信号入射角度
w = [pi/6,pi/4];%角频率
w = w';
snapshots = 200;%快拍数
D = length(w);%信号源数目
M = 12;%天线阵元数目
SNR = 10;%信噪比为10dB
A1 = zeros(D,M);%构建D行M列矩阵
A2 = zeros(D,M);
for k = 1:D
      A1(k,:) = exp(-1i*2*pi*d1*sin(theta(k))/lamda*[0:M-1]);
      A2(k,:) = exp(-1i*2*pi*d2*sin(theta(k))/lamda*[0:M-1]); 
end

A1 = A1.';%获得方向矩阵
A2 = A2.';

S = 4*exp(1i*(w*[1:snapshots]));%获得仿真信号(根据窄带信号的定义式)
X1 = awgn(A1*S,SNR,'measured');%接收信号
X2 = awgn(A2*S,SNR,'measured');%接收信号
Rx1 = X1*X1'/snapshots;%接收信号的协方差矩阵
Rx2 = X2*X2'/snapshots;
[Ve1,Va1] = eig(Rx1);%获取协方差矩阵的特征值Va以及特征向量Ve
[Ve2,Va2] = eig(Rx2);
En1 = Ve1(:,1:M-D);%取出M-D个对应特征值为噪声方差的特征向量
En2 = Ve2(:,1:M-D);

theta1 = -90:0.5:90;


%使用两个for循环进行嵌套，将各个角度依次带入进行谱峰搜索
for a = 1:length(theta1)
    AA1 = zeros(1,M);
    AA2 = zeros(1,M);
    for b = 0:M-1
        AA1(:,b+1) = exp(-1*1i*2*pi*d1*sin(theta1(a)/180*pi)/lamda*b); %注：在matlab中三角函数计算使用地是弧度值
        AA2(:,b+1) = exp(-1*1i*2*pi*d2*sin(theta1(a)/180*pi)/lamda*b);
    end
   AA1 = AA1.';
   AA2 = AA2.';
    P1 = AA1'*En1*En1'*AA1;
    P2 = AA2'*En2*En2'*AA2;
    Pmusic1(a) = abs(1/P1);%谱函数
    Pmusic2(a) = abs(1/P2);
end

Pmusic1 = 10*log10(Pmusic1/max(Pmusic1));%将获得的谱函数归一化后进行计算dB操作
Pmusic2 = 10*log10(Pmusic2/max(Pmusic2));
plot(theta1,Pmusic1,'-r')
hold on
plot(theta1,Pmusic2,'.-k')
hold off
xlabel('角度 \theta/degree')
ylabel('谱函数 P(\theta)/dB')
title('基于MUSIC算法的DOA估计(利用两个互质阵列消除由于阵元间距大于半波长引起的角度模糊)')
legend('d1=m1*lamda/2,(m1=5)','d2=m2*lamda/2,(m2=6)')
grid on

