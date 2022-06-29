%% 频响曲线分析
clc; clear all; close all;

%% 简单的例子
% y[n] = 0.9*y[n-1] + 0.1*x[n]
% h = 0.1 / (1 - 0.9*z^-1)
b = [0.1 0];
a = [1 -0.9];
[h, w] = freqz(b, a, 'whole', 2001);

figure;
subplot(2,1,1)
plot(w/pi,abs(h),'.-')
axis([0 1 -1 2])
legend('Response','Passband','Location','South')
ylabel('Magnitude')

subplot(2,1,2)
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')


%% 基于衰减曲线的混响估算
nFreq = 256;
fs = 16000;
lenc = 10;
iw = 1 : 1 : lenc;
T60 = 1.5; % 基础混响时间 
T60db = T60 * ones(nFreq, 1);         % 追踪晚期混响PSD的变量
coeff = (1 - 0.5/nFreq * (1:nFreq)).'; % 令基础混响时间随频率升高而升高
T60db = coeff .* T60db; 
alpha = 3*log(10)./(T60db*fs);  
W = exp(-2 * alpha * nFreq.*iw) / lenc;

ws = W(1, :);
b = [1 0 -1];
a = [1 0];
[h, w] = freqz(b, a, 'whole', 2001);

figure;
subplot(2,1,1)
plot(w/pi,abs(h),'.-')
axis([0 1 -1 2])
legend('Response','Passband','Location','South')
ylabel('Magnitude')

subplot(2,1,2)
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

%% 基于衰减曲线的混响估算
n=100;
w=hanning(n);
fvtool(w);
%%
b = [1 -1];
a = [1 0];
[h, w] = freqz(b, a, 'whole', 2001);

figure;
subplot(2,1,1)
plot(w/pi,abs(h),'.-')
axis([0 1 -1 2])
legend('Response','Passband','Location','South')
ylabel('Magnitude')

subplot(2,1,2)
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

%%
b = [1 -0.5];
a = [-0.5 1];
[h, w] = freqz(b, a, 'whole', 2001);

figure;
subplot(2,1,1)
plot(w/pi,abs(h),'.-')
axis([0 1 -1 2])
legend('Response','Passband','Location','South')
ylabel('Magnitude')

subplot(2,1,2)
plot(w/pi,angle(w),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')
%%
b = [0.4142 0 -0.4142];
a = [1.4142  -1.8478 1-0.4142];
[h, w] = freqz(b, a, 'whole', 2001);

figure;
subplot(2,1,1)
plot(w/pi,abs(h),'.-')
axis([0 1 -1 2])
legend('Response','Passband','Location','South')
ylabel('Magnitude')

subplot(2,1,2)
plot(w/pi,angle(w),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')
%%
b = [1 2 1];
a = [1 0 0.1716];
[h, w] = freqz(b, a, 'whole', 2001);

figure;
subplot(2,1,1)
plot(w/pi,abs(h),'.-')
axis([0 1 -1 2])
legend('Response','Passband','Location','South')
ylabel('Magnitude')

subplot(2,1,2)
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')