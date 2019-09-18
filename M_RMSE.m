%阵元个数与RMSE的关系
clc
clear
snr=10;%信噪比(dB)
sensor_number=2:12;%阵元数目
N1 = length(sensor_number);
N2 = 300;
bbb=zeros(1,N1);%用于存储各个信噪比下的均方根误差
aaa=zeros(1,N2);%用于存储估计的波达方向角度

for kk=1:N1
    
    
    for k=1:N2
        
        source_number=1;%信号源数目
        
        N_x=1024;%快拍数
        w=pi/4;%角频率
        l=2*pi*3e8/w;%波长lamda
        d=0.5*l;%阵元间距
        sigma = 0.4;
        
       
        
        source_doa=50;%信号源的波达方向
        A=[exp(-1i*(0:sensor_number(kk)-1)*d*2*pi*sin(source_doa*pi/180)/l)];%构造方向矩阵
        A = A.';%作转置运算
        
        
        s=sqrt(sigma*10.^(snr/10))*exp(1i*w*[0:N_x-1]);
        %         x = awgn(A*s,snr,'measured');
        x = A*s + (randn(sensor_number(kk),N_x)+1i*randn(sensor_number(kk),N_x));
        
        R=x*x'/N_x;
         
        [V,D]=eig(R);%对X的协方差矩阵进行特征值分解（V为特征向量，D为特征值）
        D=diag(D);
        disp(D);
        Un=V(:,1:sensor_number(kk)-source_number);%噪声部分所对应的特征向量构成的噪声子空间
        Gn=Un*Un';
        
        searching_doa=-90:0.1:90;%
        
        for i=1:length(searching_doa)
            a_theta=exp(-1i*(0:sensor_number(kk)-1)'*2*pi*d*sin(pi*searching_doa(i)/180)/l);%构造阵列流型
            Pmusic(i)=1./abs((a_theta)'*Gn*a_theta);%构造谱函数
        end
        
        %找出使得Pmusic最大时所对应的角度值
        %----找出Pmusic这组数中的极大值-----
        aa=diff(Pmusic);%Pmusic事实上是1*length(searching_doa)矩阵，用后一个数依次减去前一个数得到aa
        aa=sign(aa);%获得符号函数（aa<0时变为-1，aa=0时变为0，aa>0时取为1）
        aa=diff(aa);
        bb=find(aa==-2)+1;
        %-------------------------------------------
        
        [t1,t2]=max(Pmusic(bb));
        estimated_source_doa=searching_doa(bb(t2));
        aaa(:,k)=estimated_source_doa;
        
        %倘若只有一个信号源，可以直接使用下列代码段，但是在实际应用中，信号源数并不是已知的
        %          [a1,a2]=max(Pmusic);
        %         estimated_source_doa=searching_doa(a2);
        %         aaa(:,k)=estimated_source_doa;
        
    end
    disp(aaa);
    
  
    E_source_doa=sum(aaa(1,:))/N2;
    disp('E_source_doa');
    disp(E_source_doa);
    
    RMSE_source_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/N2);%均方根误差
    disp('RMSE_source_doa');
    disp(RMSE_source_doa);
    
    bbb(:,kk)=RMSE_source_doa;
end
disp(bbb);

plot(sensor_number,bbb(1,:),'k*-');
hold on

%CRB克拉美罗界（一维DOA估计）---阵元数M

clc
clear

snr=10;%信噪比(dB)
source_number=1;%信号源数目
sensor_number=2:12;%阵元数目
N_x=1024;%快拍数
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
hold off
xlabel('阵元数M')
ylabel('RMSE(CRB)')
title('RMSE(CRB)与阵元数关系曲线')
legend('RMSE','CRB')
grid on












