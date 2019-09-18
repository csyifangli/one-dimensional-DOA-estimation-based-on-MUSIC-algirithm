clear all;
clc;


%CLASSIC MUSIC RMSE
bbb=zeros(1,10);
for kk=1:10
sensor_number=[3 4 6 8 10 12 14 16 18 20];%以阵元数目为自变量

aaa=zeros(1,300); %每一个阵元数下，都做300次MonteCarlo试验
for k=1:300

N_x=1024;
snapshot_number=N_x;%快拍数
w=pi/4;
l=2*pi*3e8/w;
d=0.5*l;%阵元间距
snr=0;%信噪比(dB)
source_number=1; %一个信号源
source_doa=50;%信号源的入射角度
A=[exp(-j*(0:sensor_number(kk)-1)*d*2*pi*sin(source_doa*pi/180)/l)].';
s=10.^((snr/2)/10)*exp(j*w*[0:N_x-1]);
x=A*s+(1/sqrt(2))*(randn(sensor_number(kk),N_x)+j*randn(sensor_number(kk),N_x));%加了高斯白噪声后的接收数据矢量
R=x*x'/snapshot_number;


[V,D]=eig(R);
D=diag(D);

Un=V(:,1:sensor_number(kk)-source_number); %噪声子空间
Gn=Un*Un';

searching_doa=-90:0.1:90;
 for i=1:length(searching_doa)
   a_theta=exp(-j*(0:sensor_number(kk)-1)'*2*pi*d*sin(pi*searching_doa(i)/180)/l);
   Pmusic(i)=1./abs((a_theta)'*Gn*a_theta); %构建经典MUSIC的空间谱
 end

%找出谱峰所对应的极值点，即估计出入射角度
aa=diff(Pmusic);
aa=sign(aa);
aa=diff(aa);
bb=find(aa==-2)+1;
[t1 t2]=max(Pmusic(bb));
estimated_source_doa=searching_doa(bb(t2));

aaa(:,k)=estimated_source_doa;

end  %内层循环结束

%求解均方根误差
RMSE_source_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/300);%在该阵元数下，做300次试验的均方根误差
bbb(:,kk)=RMSE_source_doa;

end  %外层循环结束


plot(sensor_number,bbb(1,:),'rs-');  %绘制MUSIC算法RMSE随阵元数变化的曲线
save CLASSICAL_MUSIC_zhenyuan.mat;




%ROOT MUSIC RMSE
bbb=zeros(1,10);
for kk=1:10
sensor_number=[3 4 6 8 10 12 14 16 18 20];  %以阵元数目为自变量

aaa=zeros(1,300);
for k=1:300
    
A=[exp(-j*(0:sensor_number(kk)-1)*d*2*pi*sin(source_doa*pi/180)/l)].';
s=10.^((snr/2)/10)*exp(j*w*[0:N_x-1]);
x=A*s+(1/sqrt(2))*(randn(sensor_number(kk),N_x)+j*randn(sensor_number(kk),N_x));
R=(x*x')/snapshot_number;

[V,D]=eig(R);%对协方差矩阵进行特征分解
Un=V(:,1:sensor_number(kk)-source_number);
Gn=Un*Un'
%找出多项式的系数，并按阶数从高至低排列
a = zeros(2*sensor_number(kk)-1,1)';
for i=-(sensor_number(kk)-1):(sensor_number(kk)-1)
    a(i+sensor_number(kk)) = sum( diag(Gn,i) );
end

%使用ROOTS函数求出多项式的根
a1=roots(a)
%找出在单位圆里且最接近单位圆的N个根
a2=a1(abs(a1)<1);
%挑选出最接近单位圆的N个根
[lamda,I]=sort(abs(abs(a2)-1));
f=a2(I(1:source_number));
%计算信号到达方向角
estimated_source_doa=[asin(angle(f(1))/pi)*180/pi];
aaa(:,k)=estimated_source_doa;

end  %内层循环结束

%求解均方根误差
RMSE_source_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/300);%在该阵元数下，做300次试验的均方根误差
bbb(:,kk)=RMSE_source_doa;

end %外层循环结束


hold on
plot(sensor_number,bbb(1,:),'gp-'); %在同一幅图中绘制ROOT-MUSIC算法RMSE随阵元数变化的曲线

save ROOT_MUSIC_zhenyuan.mat;

legend('经典MUSIC','ROOT-MUSIC');
xlabel('阵元数目/个');
ylabel('估计均方根误差/度');
grid on;