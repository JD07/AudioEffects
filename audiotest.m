%% AudioDelay测试
clc; clear all; close all;

%% 
fs = 16000;
fc = 4000;
Q = 10;
audio = AudioFilter();
audio.reset(fs);
params = audio.getParameters();
params.algorithm = "kResonB";
params.fc = fc;
params.Q = Q;
params.boostCut_dB = 10;
audio.setParameters(params);

a0 = audio.coeffArray(1);
a1 = audio.coeffArray(2);
a2 = audio.coeffArray(3);
b1 = audio.coeffArray(4);
b2 = audio.coeffArray(5);

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