%Root-MUSIC算法的RMSE曲线

clear
clc

bbb=zeros(1,8);
for kk=1:8
    snr=[-5 -3 0 2 4 6 8 10]
    aaa=zeros(2,300);
    for k=1:300             %Monte carno试验次数
        Monte=300;
        kelm=16;
        derad = pi/180;
        radeg = 180/pi;
        twpi = 2*pi;
        dd = 0.5;
        d=0:dd:(kelm-1)*dd;
        iwave = 2;              % number of DOA
        theta = [10 30];
        n =200;
        A=exp(-1i*twpi*d.'*(sin(theta*derad)));
        S=randn(iwave,n);
        X0=A*S;
        X=awgn(X0,snr(kk),'measured');
        Rxx=X*X';
        InvS=inv(Rxx);
        [EVx,Dx]=eig(Rxx);
        EVAx=diag(Dx)';
        [EVAx,Ix]=sort(EVAx);
        EVAx=fliplr(EVAx);
        EVx=fliplr(EVx(:,Ix));
        Unx=EVx(:,iwave+1:kelm);
        syms z
        pz = z.^([0:kelm-1]');
        pz1 = (z^(-1)).^([0:kelm-1]);
        fz = z.^(kelm-1)*pz1*Unx*Unx'*pz;
        a = sym2poly(fz)
        zx = roots(a)
        rx=zx.';
        [as,ad]=(sort(abs((abs(rx)-1))));
        DOAest=asin(sort(-angle(rx(ad([1,3])))/pi))*180/pi;
        aaa(:,k)=DOAest.';
    end
    
    
    RMSE1=sqrt(sum((aaa(1,:)-theta(1,1)).^2)/Monte);
    RMSE2=sqrt(sum((aaa(2,:)-theta(1,2)).^2)/Monte);
    RMSE=(RMSE1+RMSE2)/2;
    bbb(:,kk)=RMSE;
end

plot(snr,bbb(1,:),'k*-');
xlabel('信噪比SNR/dB')
ylabel('均方根误差RMSE/degree')
title('经典MUSIC算法DOA估计的RMSE曲线')
grid on


