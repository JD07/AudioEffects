%% 一阶二阶滤波器仿真
% 参考<Designing Audio Effect Plugins in C++>
clc; clear all; close all;

%% Simple Resonator
% 谐振器就是一种带通滤波器，只不过其峰值十分狭窄
% 差分方程为：y[n] = a0*x[n] - b1*y[n-1] - b2*y[n-2]
% 系数计算参考（11.4）
fc = 2000;
fs = 48000;
Q = 5;
theta = 2*pi*fc/fs;
BW = fc/Q;
b2 = exp(-2*pi*BW/fs);
b1 = -4*b2*cos(theta)/(1+b2);
a0 = (1-b2)*sqrt(1-b1^2/(4*b2));
a1 = 0; a2 = 0;

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

%% Smith-Angell Resonator
% 该谐振器在z=1和z=-1处增加零点，令滤波器更加陡峭
% y[n] = a0*x[n] + a2*x[n-2] - b1*y[n-1] -b2*y[n-2]
fc = 2000;
fs = 48000;
theta = 2*pi*fc/fs;
BW = fc/Q;
b2 = exp(-2*pi*BW/fs);
b1 = -4*b2*cos(theta)/(1+b2);
a0 = 1-sqrt(b2);
a2 = -a0;
a1 = 0;

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

%% 1st Order LPF and HPF (11.24)
fc = 2000;
fs = 16000;
theta = 2*pi*fc / fs;
gamma = cos(theta)/(1 + sin(theta));

% LPF
a0 = (1-gamma)/2;
a1 = (1-gamma)/2;
a2 = 0.0;
b1 = -gamma;
b2 = 0.0;

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,abs(h),'.-')
axis([0 1 0 1])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

% HPF
a0 = (1+gamma)/2;
a1 = -(1+gamma)/2;
a2 = 0.0;
b1 = -gamma;
b2 = 0.0;

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,abs(h),'.-')
axis([0 1 0 1])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

%% 2nd Order LPF and HPF (11.25)
fc = 2000;
fs = 16000;
Q = 0.707;
theta = 2*pi*fc/fs;
d = 1/Q;
beta = 0.5 * (1 - 0.5*d*sin(theta))/(1 + 0.5*d*sin(theta));
gamma = (0.5 + beta)*cos(theta);

% LPF
a0 = (0.5 + beta - gamma)/2;
a1 = 0.5 + beta - gamma;
a2 = (0.5 + beta - gamma)/2;
b1 = -2*gamma;
b2 = 2*beta;

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,abs(h),'.-')
axis([0 1 0 1])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

% HPF
a0 = (0.5 + beta + gamma)/2;
a1 = -(0.5 + beta + gamma);
a2 = (0.5 + beta + gamma)/2;
b1 = -2*gamma;
b2 = 2*beta;

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,abs(h),'.-')
axis([0 1 0 1])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

%% 2nd Order BPF and BSF (11.26)
fc = 4000;
fs = 16000;
Q = 0.707;
K = tan(pi*fc/fs);
sigma = K^2 * Q + K + Q;

% BPF
a0 = K/sigma;
a1 = 0.0;
a2 = -K/sigma;
b1 = 2*Q*(K^2 - 1)/sigma;
b2 = (K^2*Q - K + Q)/sigma;

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,abs(h),'.-')
axis([0 1 0 1])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

% BSF
a0 = Q*(K^2 + 1)/sigma;
a1 = 2*Q*(K^2 - 1)/sigma;
a2 = Q*(K^2 + 1)/sigma;
b1 = 2*Q*(K^2 - 1)/sigma;
b2 = (K^2 * Q - K + Q)/sigma;
a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,abs(h),'.-')
axis([0 1 0 1])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

%% 2nd Order Butterworth LPF and HPF (11.27)
% Butterworth是特殊的2阶LPF和HPF，其相当于Q被固定于0.707
% Q=0.707是观察不到突出峰的极限值
fc = 4000;
fs = 16000;

% LPF
C = 1/tan(pi*fc/fs);
a0 = 1/(1 + sqrt(2)*C + C^2);
a1 = 2*a0;
a2 = a0;
b1 = 2*a0*(1 - C^2);
b2 = a0*(1 - sqrt(2)*C + C^2);

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,abs(h),'.-')
axis([0 1 0 1])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

% HPF
C = 1/tan(pi*fc/fs);
a0 = 1/(1 + sqrt(2)*C + C^2);
a1 = -2*a0;
a2 = a0;
b1 = 2*a0*(C^2 - 1);
b2 = a0*(1 - sqrt(2)*C + C^2);

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,abs(h),'.-')
axis([0 1 0 1])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

%% 2nd Order Butterworth BPF and BSF (11.28)
fc = 4000;
fs = 16000;
Q = 5;
BW = fc/Q;

% BPF
C = 1/tan(pi*BW/fs);
D = 2*cos(2*pi*fc/fs);
a0 = 1/(1+C);
a1 = 0.0;
a2 = -a0;
b1 = -a0*C*D;
b2 = a0*(C - 1);

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,abs(h),'.-')
axis([0 1 0 1])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

% BSF
C = tan(pi*BW/fs);
D = 2*cos(2*pi*fc/fs);
a0 = 1/(1+C);
a1 = -a0*D;
a2 = a0;
b1 = -a0*D;
b2 = a0*(1 - C);

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,abs(h),'.-')
axis([0 1 0 1])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

%% 2nd Order Linkwitz-Riley LPF and HPF (11.29)
% 2阶Linkwitz-Riley滤波器在拐点处的衰减频率为-6dB而不是-3dB
fc = 4000;
fs = 16000;
theta = pi*fc/fs;
omega = pi*fc;
k = omega/tan(theta);
sigma = k^2 + omega^2 + 2*k*omega;

% LPF
a0 = omega^2/sigma;
a1 = 2*omega^2/sigma;
a2 = omega^2/sigma;
b1 = (-2*k^2 + 2*omega^2)/sigma;
b2 = (-2*k*omega + k^2 + omega^2)/sigma;

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,abs(h),'.-')
axis([0 1 0 1])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

% HPF
a0 = k^2/sigma;
a1 = -2*k^2/sigma;
a2 = k^2/sigma;
b1 = (-2*k^2 + 2*omega^2)/sigma;
b2 = (-2*k*omega + k^2 + omega^2)/sigma;

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,abs(h),'.-')
axis([0 1 0 1])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

%% 1st and 2nd Order APF (11.30)
fc = 4000;
fs = 16000;

% 1st Order APF
alpha = (tan(pi*fc/fs) - 1)/(tan(pi*fc/fs) + 1);
a0 = alpha;
a1 = 1.0;
a2 = 0.0;
b1 = alpha;
b2 = 0.0;

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,abs(h),'.-')
axis([0 1 0 1])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

% 2st Order APF
Q = 2;
BW = fc / Q;
alpha = (tan(pi*BW/fs) - 1)/(tan(pi*BW/fs) + 1);
beta = -cos(2*pi*fc/fs);
a0 = -alpha;
a1 = beta*(1-alpha);
a2 = 1.0;
b1 = beta*(1-alpha);
b2 = -alpha;

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,abs(h),'.-')
axis([0 1 0 1])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

%% 1st Order Low Shelving and High Shelving
% 注意，这个low-shelving的实现和另一个教材不一样
% 其本质上实现了一个低通/高通滤波器，然后滤波器输出乘以倍数后和原信号进行叠加
fc = 4000;
fs = 16000;
G_dB = 10; %增益，单位dB
theta = 2*pi*fc/fs;
u = 10^(G_dB/20);

% low shelving
beta = 4/(1 + u);
sigma = beta*tan(theta/2);
gamma = (1-sigma)/(1+sigma);
a0 = (1-gamma)/2;
a1 = (1-gamma)/2;
a2 = 0.0;
b1 = -gamma;
b2 = 0.0;

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,(u-1)*abs(h) + 1,'.-')
axis([0 1 0 u])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

% high shelving
beta = (1 + u)/4;
sigma = beta*tan(theta/2);
gamma = (1-sigma)/(1+sigma);
a0 = (1+gamma)/2;
a1 = -(1+gamma)/2;
a2 = 0.0;
b1 = -gamma;
b2 = 0.0;

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,(u-1)*abs(h) + 1,'.-')
axis([0 1 0 u])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

%% 2nd Order Parametric EQ Filter: Non-Constant-Q
fc = 4000;
fs = 16000;
G_dB = 10; %增益，单位dB
theta = 2*pi*fc/fs;
u = 10^(G_dB/20);
sigma = 4/(1 + u);
beta = 0.5 * (1 - sigma*tan(theta/(2*Q))) / (1 + sigma*tan(theta/(2*Q)));
gamma = (0.5 + beta)*cos(theta);
a0 = 0.5-beta;
a1 = 0.0;
a2 = -(0.5-beta);
b1 = -2*gamma;
b2 = 2*beta;

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,(u-1)*abs(h) + 1,'.-')
axis([0 1 0 u])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

%% 2nd Order Parametric EQ Filter:Constant-Q
fc = 4000;
fs = 16000;
Q = 1;
G_dB = 10; %增益，单位dB
K = tan(pi*fc/fs);
Vo = 10^(G_dB/20);
d0 = 1 + K/Q + K^2;
e0 = 1 + K/(Vo*Q) + K^2;

alpha = 1 + Vo*K/Q + K^2;
beta = 2*(K^2 - 1);
gamma = 1 - Vo*K/Q + K^2;
delta = 1 - K/Q + K^2;
eta = 1 - K/(Vo*Q) + K^2;

% Boost
a0 = alpha/d0;
a1 = beta/d0;
a2 = gamma/d0;
b1 = beta/d0;
b2 = delta/d0;

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,(u-1)*abs(h) + 1,'.-')
axis([0 1 0 u])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

% Cut
a0 = d0/e0;
a1 = beta/e0;
a2 = delta/e0;
b1 = beta/e0;
b2 = eta/e0;

a = [a0 a1 a2];
b = [1 b1 b2];
[h, w] = freqz(a, b, 'whole', 2001);
figure();
subplot(2,1,1);
plot(w/pi,(u-1)*abs(h) + 1,'.-')
axis([0 1 0 u])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Gain')
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')
