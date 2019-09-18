clc
clear

Ppairs = [0.0045448,0.0151145,0.0188747,0.0390625];
Rpairs = [0.0165566,0.0235808,0.0265344,0.0312500];



for i = 1:4
    F1 = 2*Ppairs(i)*Rpairs(i)/(Rpairs(i)+Ppairs(i));
    disp(F1)
end

