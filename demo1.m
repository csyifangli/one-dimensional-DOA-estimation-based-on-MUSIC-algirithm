clc
clear all
format long %the data show that as long shaping scientific
doa = [20 60]/180*pi;%direction of arrival
N = 200;%snapshots
w = [pi/4 pi/3]';%Frequency
M = 10;%number of array elements
P = length(w);%the number of signal
lambda = 150;%wavelength
d = lambda/2;%element spacing
snr = 20;%SNR
D = zeros(P,M);%create a matrix with P rows and M column
for k = 1:P
    D(k,:) = exp(-j*2*pi*d*sin(doa(k))/lambda*[0:M-1]);%assignment matrix
end
D = D';

% E = zeros(M,P);
% for k = 1:P
%     E(:,k) = exp(-j*2*pi*d*sin(doa(k))/lambda*[0:M-1]);%assignment matrix
% end
% 
%     disp(D);
%     disp(E);


xx = 2*exp(j*(w*[1:N]));%simulate signal(narrowband signal)
x = D*xx;
x = x + awgn(x,snr);%Insert Gaussian white noise
R  = x*x';%data covariance matrix
[N,V] = eig(R);%find the eigenvalues and eigenvectors of R
NN = N(:,1:M-P);%estimate noise subspace
theta = -90:0.5:90;%peak search
for ii = 1:length(theta)
    SS = zeros(1,M);
    for jj = 0:M-1
        SS(1+jj) = exp(-j*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda);
    end
    PP = SS*NN*NN'*SS';%get a concrete number
    Pmusic(ii) = abs(1/PP);
end
Pmusic = 10*log10(Pmusic/max(Pmusic));%spatial spectrum function
plot(theta,Pmusic,'-k')
xlabel('angle \theta/degree')
ylabel('spectrum function P(\theta) /dB')
title('DOA estimation based on MUSIC algorithm')
grid on




