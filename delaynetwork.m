% delaynetwork. m
% Author: T. Lokki
fs=44100;
gain=0.99;
% Create an impulse
x = zeros(1, 1*fs); x(1) = 1;
y = zeros(1, fs);
b = [1 1 1 1];
c = [0.8 0.8 0.8 0.8];
% Feedback matrix
% 反馈矩阵是一个酉矩阵
a(1, : ) =[0 1 1 0];
a(2, : ) =[-1 0 0 -1];
a(3, : ) =[1 0 0 -1];
a(4, : ) =[0 1 -1 0];
a2=a*(1/sqrt(2) ) * gain;
% Delay lines, use prime numbers
% 使用素数，防止不同延时之间产生相干性，导致结果染色
m=[149 211 263 293]';
z1=zeros(1, max(max(m) ) );
z2=zeros(1, max(max(m) ) );
z3=zeros(1, max(max(m) ) );
z4=zeros(1, max(max(m) ) );
for n = 1: length(y)
tmp = [z1(m(1) ) z2(m(2) ) z3(m(3) ) z4(m(4) ) ];
y(n) = x(n) + c(1) *z1(m(1) ) + c(2) *z2(m(2) ) ...
+ c(3) *z3(m(3) ) + c(4) *z4(m(4) );
z1 = [(x(n)*b(1) + tmp*a2(1, : )') z1(1: length(z1) -1)];
z2 = [(x(n)*b(2) + tmp*a2(2, : )') z2(1: length(z2) -1)];
z3 = [(x(n)*b(3) + tmp*a2(3, : )') z3(1: length(z3) -1)];
z4 = [(x(n)*b(4) + tmp*a2(4, : )') z4(1: length(z4) -1)];
end
figure
plot(y)

%% 
