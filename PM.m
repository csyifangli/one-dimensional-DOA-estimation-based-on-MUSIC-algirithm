%利用PM传播因子获取噪声子空间，避免特征值分解
%利用代价函数去估计传播算子P
clc
clear

lamda = 100;%波长设为100
d = lamda/2;%阵元间距取为半波长
theta = [20,30,60];%信号源入射角度
M = 12;%阵元个数
SNR = 10;%信噪比为10dB
w = [pi/6,pi/4,pi/3];%角频率
w = w';
snapshots = 1024;%快拍数
D = length(theta);%信号源数

%构造方向矩阵
A = zeros(M,D);

for k = 1:D
    A(:,k) = exp(-1i*2*pi*d*sin(theta(k)*pi/180)*[0:M-1]/lamda);
end

S = 4*exp(1i*(w*[1:snapshots]));%获得仿真信号(根据窄带信号的定义式)
X = awgn(A*S,SNR,'measured');%接收信号
Rx = X*X'/snapshots;%接收信号的协方差矩阵

%获取传播算子
%将协方差矩阵进行分块
G = Rx(:,1:D);
H = Rx(:,D+1:M);

%根据在有噪声情况下通过代价函数获取的传播算子的估计值

P = inv(G'*G)*G'*H;

Q = [P;-eye(M-D)];%噪声子空间

%使用其正交化后的Q来取代原来的Q，使其列向量正交，更加逼近原来的噪声子空间
Q0 = Q*(Q'*Q)^(-1/2);

Gn = Q0*Q0';


Stheta = -90:0.5:90;
    
for k1 = 1:length(Stheta)
    AA = zeros(1,M);
    for k2 = 0:M-1
        AA(:,k2+1) = exp(-1i*2*pi*d*sin(Stheta(k1)*pi/180)/lamda*k2);
    end
    
    AA = AA.';
    Pmusic(k1) = 1/abs(AA'*Gn*AA);
    
end

%画图
Pmusic = 10*log10(Pmusic/max(Pmusic));
plot(Stheta,Pmusic,'-r')
hold on
xlabel('角度 \theta/degree')
ylabel('谱函数 P(\theta)/dB')
title('基于MUSIC算法的DOA估计(PM替代特征值分解获取噪声子空间)')
axis([-90,90,-50,0])
grid on

