%the relationship between DOA estimation and the number of array elements

clc
clear all
format long %The data show that as long shaping science and technology
N=200; %Snapshots
doa=[20 60]/180*pi; %DOA
w=[pi/4 pi/3]'; %Frequency
M1=10; %Array elements number
M2=50;
M3=100;
P=length(w); %Number of signal
lambda=150; %Wavelength
d=lambda/2; %Array element spacing
snr=20; %SNR
D1=zeros(P,M1);
D2=zeros(P,M2);
D3=zeros(P,M3);
for k=1:P
D1(k,:)=exp(-j*2*pi*d*sin(doa(k))/lambda*[0:M1-1]); %Assignment matrix
D2(k,:)=exp(-j*2*pi*d*sin(doa(k))/lambda*[0:M2-1]);
D3(k,:)=exp(-j*2*pi*d*sin(doa(k))/lambda*[0:M3-1]);
end
D1=D1';
D2=D2';
D3=D3';
xx=2*exp(j*(w*[1:N])); %Simulate the signal
x1=D1*xx;
x2=D2*xx;
x3=D3*xx;
x1=x1+awgn(x1,snr); %Add Gaussian white noise
x2=x2+awgn(x2,snr);
x3=x3+awgn(x3,snr);
R1=x1*x1'; %Data covariance matrix
R2=x2*x2';
R3=x3*x3';
[N1,V1]=eig(R1); %Find the eigenvalues and eigenvectors of R
[N2,V2]=eig(R2);
[N3,V3]=eig(R3);
NN1=N1(:,1:M1-P); ; %Estimate the noise subspace
NN2=N2(:,1:M2-P);
NN3=N3(:,1:M3-P);
theta=-90:0.5:90;
%%Search the peak
for ii=1:length(theta)
SS1=zeros(1,M1);
for jj=0:M1-1
SS1(1+jj)=exp(-j*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda);
end
PP1=SS1*NN1*NN1'*SS1';
Pmusic1(ii)=abs(1/ PP1);
end
for ii=1:length(theta)
SS2=zeros(1,M2);
for jj=0:M2-1
SS2(1+jj)=exp(-j*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda);
end
PP2=SS2*NN2*NN2'*SS2';
Pmusic2(ii)=abs(1/ PP2);
end
for ii=1:length(theta)
SS3=zeros(1,M3);
for jj=0:M3-1
SS3(1+jj)=exp(-j*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda);
end
PP3=SS3*NN3*NN3'*SS3';
Pmusic3(ii)=abs(1/ PP3);
end
Pmusic1=10*log10(Pmusic1/max(Pmusic1)); %Spatial spectrum function
Pmusic2=10*log10(Pmusic2/max(Pmusic2));
Pmusic3=10*log10(Pmusic3/max(Pmusic3));
plot(theta,Pmusic1,'--k','LineWidth',2.0)
hold on
plot(theta,Pmusic2,'k','LineWidth',1.0)
hold on
plot(theta,Pmusic3,'-.k','LineWidth',0.1)
hold off
xlabel('angle \theta/degree')
ylabel('spectrum function P(\theta) /dB')
title('DOA estimation based on MUSIC algorithm')
legend('10 array elements','50 array elements','100 array elements')
grid on