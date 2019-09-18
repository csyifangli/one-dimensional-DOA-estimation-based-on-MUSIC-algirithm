%������snapshots��RMSE�Ĺ�ϵ
clc
clear
snr=0;%�����(dB)
N_x=100:200:4096;%������
N1 = length(N_x);
N2 = 300;
bbb=zeros(1,N1);%���ڴ洢����������µľ��������
aaa=zeros(1,N2);%���ڴ洢���ƵĲ��﷽��Ƕ�

for kk=1:N1
    
    
    for k=1:N2
        
        source_number=1;%�ź�Դ��Ŀ
        sensor_number=8;%��Ԫ��Ŀ
        
        w=pi/4;%��Ƶ��
        l=2*pi*3e8/w;%����lamda
        d=0.5*l;%��Ԫ���
        sigma = 0.4;
        
        
        source_doa=50;%�ź�Դ�Ĳ��﷽��
        A=[exp(-1i*(0:sensor_number-1)*d*2*pi*sin(source_doa*pi/180)/l)];%���췽�����
        A = A.';%��ת������
        
      
        s=sqrt(sigma*10.^(snr/10))*exp(1i*w*[0:N_x(kk)-1]);
%         x = awgn(A*s,snr,'measured');
        x = A*s + (randn(sensor_number,N_x(kk))+1i*randn(sensor_number,N_x(kk)));
        
        R=x*x'/N_x(kk);
        
        [V,D]=eig(R);%��X��Э��������������ֵ�ֽ⣨VΪ����������DΪ����ֵ��
        D=diag(D);%ȡ��D�Խ���Ԫ��
        disp(D);
        Un=V(:,1:sensor_number-source_number);%������������Ӧ�������������ɵ������ӿռ�
        Gn=Un*Un';
        
        searching_doa=-90:0.1:90;%
        
        for i=1:length(searching_doa)
            a_theta=exp(-1i*(0:sensor_number-1)'*2*pi*d*sin(pi*searching_doa(i)/180)/l);%������������
            Pmusic(i)=1./abs((a_theta)'*Gn*a_theta);%�����׺���
        end
        
        %�ҳ�ʹ��Pmusic���ʱ����Ӧ�ĽǶ�ֵ
        %----�ҳ�Pmusic�������еļ���ֵ-----
        aa=diff(Pmusic);%Pmusic��ʵ����1*length(searching_doa)�����ú�һ�������μ�ȥǰһ�����õ�aa
        aa=sign(aa);%��÷��ź�����aa<0ʱ��Ϊ-1��aa=0ʱ��Ϊ0��aa>0ʱȡΪ1��
        aa=diff(aa);
        bb=find(aa==-2)+1;
        %-------------------------------------------
        
        [t1,t2]=max(Pmusic(bb));
        estimated_source_doa=searching_doa(bb(t2));
        aaa(:,k)=estimated_source_doa;
        
        %����ֻ��һ���ź�Դ������ֱ��ʹ�����д���Σ�������ʵ��Ӧ���У��ź�Դ����������֪��
        %          [a1,a2]=max(Pmusic);
        %         estimated_source_doa=searching_doa(a2);
        %         aaa(:,k)=estimated_source_doa;
        
    end
    disp(aaa);
    
    
    E_source_doa=sum(aaa(1,:))/N2;
    disp('E_source_doa');
    disp(E_source_doa);
    
    RMSE_source_doa=sqrt(sum((aaa(1,:)-source_doa).^2)/N2);%���������
    disp('RMSE_source_doa');
    disp(RMSE_source_doa);
    
    bbb(:,kk)=RMSE_source_doa;
end
disp(bbb);

plot(N_x,bbb(1,:),'k*-');
hold on

%CRB�������޽磨һάDOA���ƣ�---������

clc
clear

snr=0;%�����(dB)
source_number=1;%�ź�Դ��Ŀ
sensor_number=8;%��Ԫ��Ŀ
N_x=100:200:4096;%������
w=pi/4;%��Ƶ��
l=2*pi*3e8/w;%����lamda
d=0.5*l;%��Ԫ���
sigma = 0.4;%��������

source_doa=[50]*pi/180;%�ź�Դ�Ĳ��﷽��
A = zeros(source_number,sensor_number);%����D��M�о���
for k = 1:source_number
      A(k,:) = exp(-1i*2*pi*d*sin(source_doa(k))/l*[0:sensor_number-1]);
end
A = A.';%��ת������

%A��theta����ΪD
D = zeros(sensor_number,source_number);
for k1 = 1:source_number
    for k2 = 0:sensor_number-1
        D(k2+1,k1) = (-1i*2*pi*d*cos(source_doa(k1))/l*k2*pi/180)*exp(-1i*2*pi*d*sin(source_doa(k1))/l*k2);
    end
end

crb = zeros(1,length(N_x));
for kk = 1:length(N_x)
    S = sqrt(sigma*10.^(snr/10))*exp(1i*w*[0:N_x(kk)-1]);
       
    P = S*S'/N_x(kk);%�����źŵ�Э�������
    
    R = A*P*A' + sigma*eye(sensor_number);
    
    crb(:,kk) = sqrt(sigma/(2*N_x(kk))*inv((real(D'*(eye(sensor_number)-A*inv(A'*A)*A')*D).*(P*A'*inv(R)*A*P).')));
    
end

%��ͼ
plot(N_x,crb(1,:), '--rd',  'LineWidth', 1.2, 'MarkerSize', 8)
hold off
xlabel('������')
ylabel('RMSE(CRB)')
title('RMSE��CRB�����������ϵ����')
grid on

