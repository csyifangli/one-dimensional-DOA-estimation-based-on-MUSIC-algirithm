%利用TLS准则去估计传播算子P
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
R = Rx'*Rx;
%对R进行特征值分解
[Ve,Va] = eig(R);

Vb = Ve(:,1:M-D);

% %对Vb进行分块处理
% V1b = Vb(1:D,:);
% V2b = Vb(D+1:M,:);
% 
% P = -V1b*inv(V2b);
% 
% Q = [P;-eye(M-D)];%噪声子空间

%此处的Q使用了上面方法的正交化来取代
Q = -Vb;
Gn = Q*Q';


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
xlabel('角度 \theta/degree')
ylabel('谱函数 P(\theta)/dB')
title('基于MUSIC算法的DOA估计(PM替代特征值分解获取噪声子空间)')
axis([-90,90,-60,0])
grid on

