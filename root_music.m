%Root-MUSIC algorithm
clc,clear
format long
lamda = 150;%wave length
d = lamda/2;%the spacing of array
L = 200;%snapshots
w = [pi/6,pi/4,pi/3];%2*pi*frequency of signal
w = w';
W = length(w);%the number of signal source
M = 10;%the number of array element
snr = -10;%SNR
theta = [20,30,60];%direction angle of signal source
theta = theta/180*pi;%angle to degree

S = 2*exp(1i*(w*[1:L]));%simulate the signal
A = zeros(M,W);
%create an assignment matrix(方向矩阵A)
for k = 1:W
    A(:,k) = exp(-1i*2*pi*d*sin(theta(k))/lamda*[0:M-1]);
end

X = awgn(A*S,snr,'measured');
R = X*X'/L;
[Ve,Va] = eig(R);

En = Ve(:,1:M-W);


%-----------------第一种多项式构造方法----P(z)=z^(M-1)*P'(z^-1)*En*En'*P(z)
G = En*En';

a = zeros(2*M-1,1);
a = a';
for k = -(M-1):(M-1) %获取多项式的系数
      a(M+k) = sum(diag(G,k));
end

a1 = roots(a);%得到多项式的解(a中存储的是降幂的多项式的系数)
%接下来从这2（M-1）个解中找出在单位圆内且最接近单位圆的W个解
a2 = a1(abs(a1)<1);%先将单位圆内的解取出
[Q,I] = sort(abs(a2),'descend');%找到离单位圆最近的解的索引I
f = a2(I(1:W));

DOA = asin(lamda/(2*pi*d)*angle(f));
DOA = sort(DOA);
disp(DOA*180/pi)

% %---------第二种方法构造的多项式---------------------------
% %
% En1 = En(1:W,:);
% En2 = En(W+1:M,:);
% B = eye(M-W);
% b = B(:,1);
% c = En1*inv(En2)*b;
% c = fliplr(c');
% 
% c1 = [1,c];
% 
% f = roots(c1);
% DOA = asin(lamda/(2*pi*d)*angle(f));
% DOA = sort(DOA);
% disp(DOA*180/pi)

