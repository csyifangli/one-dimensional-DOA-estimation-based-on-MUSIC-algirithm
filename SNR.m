%信噪比SNR与MUSIC算法的关系
clc
clear 
format long
lamda = 150;%最高频率信号的波长
d = lamda/2;%阵元间距
%此处需注意length(theta) = length(w)
theta = [20,30,60]/180*pi;%信号入射角度
w = [pi/6,pi/4,pi/3];%角频率
w = w';
snapshots = 200;%快拍数
D = length(w);%信号源数目
M = 12;%天线阵元数目
SNR = 20;%信噪比为20dB
SNR2 = 40;
SNR3 = 80;
A = zeros(D,M);%构建D行M列矩阵
for k = 1:D
      A(k,:) = exp(-1i*2*pi*d*sin(theta(k))/lamda*[0:M-1]);
end

A = A';%获得投影矩阵
S = 4*exp(1i*(w*[1:snapshots]));%获得仿真信号(根据窄带信号的定义式)
X = A*S + awgn(A*S,SNR);%接收信号
X2 = A*S + awgn(A*S,SNR2);%接收信号
X3 = A*S + awgn(A*S,SNR3);%接收信号

Rx = X*X'/snapshots;%接收信号的协方差矩阵
Rx2 = X2*X2'/snapshots;%接收信号的协方差矩阵
Rx3 = X3*X3'/snapshots;%接收信号的协方差矩阵
[Ve,Va] = eig(Rx);%获取协方差矩阵的特征值Va以及特征向量Ve
[Ve2,Va2] = eig(Rx2);
[Ve3,Va3] = eig(Rx3);



En = Ve(:,1:M-D);%取出M-D个对应特征值为噪声方差的特征向量
En2 = Ve2(:,1:M-D);
En3 = Ve3(:,1:M-D);


theta1 = -90:0.5:90;


%使用两个for循环进行嵌套，将各个角度依次带入进行谱峰搜索
for a = 1:length(theta1)
    AA = zeros(1,M);
    for b = 0:M-1
        AA(:,b+1) = exp(-1*1i*2*pi*d*sin(theta1(a)/180*pi)/lamda*b); %注：在matlab中三角函数计算使用地是弧度值
    end
    P = AA*En*En'*AA';
    Pmusic(a) = abs(1/P);%谱函数
    
    P2 = AA*En2*En2'*AA';
    Pmusic2(a) = abs(1/P2);%谱函数
    
    P3 = AA*En3*En3'*AA';
    Pmusic3(a) = abs(1/P3);%谱函数
end

Pmusic = 10*log10(Pmusic/max(Pmusic));%将获得的谱函数归一化后进行计算dB操作
Pmusic2 = 10*log10(Pmusic2/max(Pmusic2));
Pmusic3 = 10*log10(Pmusic3/max(Pmusic3));
plot(theta1,Pmusic,'-r')
hold on
plot(theta1,Pmusic2,'--k')
hoold on
plot(theta1,Pmusic3,'-.k')
hold off
xlabel('角度 \theta/degree')
ylabel('谱函数 P(\theta)/dB')
title('基于MUSIC算法的DOA估计')
legend('SNR=20','SNR=40','SNR=80')
grid on

%c = find(Pmusic == max(Pmusic));
%text(theta1(c),Pmusic(c),['(',num2str(theta1(c)),',',num2str(Pmusic(c)),')'],'color','b');






