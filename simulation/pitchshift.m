%% 模拟LFO控制延时对相位偏移的影响
clc; clear all; close all;

%% 基本参数
fs = 16000;
f = 4; % 单位：Hz
W = 0.002; % 单位：s
Mavg = 0.002; % 单位：s
n = 1:1:2*fs;

%% 正弦波
y = Mavg + W * sin(2*pi*n*f/fs);
Y1 = [0 y];
Y0 = [y 0];
fratio = (Y0 - Y1) * fs; % 注意，每个采样点对应的时间为1/fs，求导的时候不要忘了除
fratio = 1 + fratio(2:end);
%fratio = 1 - 2*pi*f*W*cos(2*pi*n*f/fs);

figure()
subplot(2,1,1);
plot(n/16000, 1000 * y);
xlabel("Time(sec)");
ylabel("Delay(ms)");

subplot(2,1,2);
semitone = 0.0595*ones(1, length(n));
plot(n/16000, fratio);
hold on;
plot(n/16000, 1 + semitone, 'r');
text(0.001, 1.08, '+1 semitone');
plot(n/16000, 1 - semitone, 'r');
text(0.001, 0.92, '-1 semitone');
axis([0 2 0.9 1.1])

xlabel("Time(sec)");
ylabel("Pitch Shift(ratio)");

%% 三角波
y =  Mavg + W * (2 * abs(2 * mod(f * n, fs)/fs -1) - 1);
Y1 = [0 y];
Y0 = [y 0];
fratio = (Y1 - Y0) * fs; % 注意，每个采样点对应的时间为1/fs，求导的时候不要忘了除
fratio = 1 + fratio(2:end);

figure()
subplot(2,1,1);
plot(n/16000, 1000 * y);
xlabel("Time(sec)");
ylabel("Delay(ms)");

subplot(2,1,2);
semitone = 0.0595*ones(1, length(n));
plot(n/16000, fratio);
hold on;
plot(n/16000, 1 + semitone, 'r');
text(0.001, 1.08, '+1 semitone');
plot(n/16000, 1 - semitone, 'r');
text(0.001, 0.92, '-1 semitone');
axis([0 2 0.9 1.1])

%% 锯齿波
y =  Mavg + W * (2 * mod(f * n, fs)/fs - 1);
Y1 = [0 y];
Y0 = [y 0];
fratio = (Y1 - Y0) * fs; % 注意，每个采样点对应的时间为1/fs，求导的时候不要忘了除
fratio = 1 + fratio(2:end);

figure()
subplot(2,1,1);
plot(n/16000, 1000 * y);
xlabel("Time(sec)");
ylabel("Delay(ms)");

subplot(2,1,2);
semitone = 0.0595*ones(1, length(n));
plot(n/16000, fratio);
hold on;
plot(n/16000, 1 + semitone, 'r');
text(0.001, 1.08, '+1 semitone');
plot(n/16000, 1 - semitone, 'r');
text(0.001, 0.92, '-1 semitone');
axis([0 2 0.9 1.1])