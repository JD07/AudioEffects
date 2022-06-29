%% 研究不同延时结构的冲激响应
clc; clear all; close all;

%% 基本参数
fs = 16000; % 采样率
delaytime = 30; % 单位ms
fb = -0.5; % 反馈系数
dry = 1.0; % 干信号系数
wet = 1.0; % 湿信号系数

%% 延时模块初始化
simpledelay = SimpleDelay();
simpledelay.reset(fs);
parameters = simpledelay.getParameters();
parameters.delayTime_mSec = delaytime;
simpledelay.setParameters(parameters);

%% 冲击响应计算
res = zeros(1, fs);
for i = 1:fs
    if 1 == i
        xn = 1;
    else
        xn = 0;
    end
    yn = simpledelay.readDelay();
    input = xn + fb * yn;
    simpledelay.processAudioSample(input);
    output = dry * xn + wet*yn;
    res(i) = output;
end

%% 结果可视化
figure();
hold on;
for i = 1:fs
    if 0 ~= res(i)
        stem(i, res(i),'or');
    end
end