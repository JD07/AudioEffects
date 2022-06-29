%% 模块测试
clc; clear; close;

%% LFO模块测试
% 2021/2/21通过测试
% 2021/2/22句柄修改为handle，通过测试

% 类的初始化以及参数设置
sampleRate = 16000;
lfo = LFO();
lfo.reset(sampleRate);
params = lfo.getParameters();
params.waveform = "kSaw";
params.frequency_Hz = 2;
lfo.setParameters(params);

% 测试结构体
sampleLength = sampleRate; % 测试时间为0.1s，因此应当有 0.1 * value2 个完整波形
normalOutput = zeros(1, sampleLength);
quadPhaseOutput_pos = zeros(1, sampleLength);
quadPhaseOutput_neg = zeros(1, sampleLength);
normalOutputoffset = zeros(1, sampleLength);
for i = 1:sampleLength
    output = lfo.renderAudioOutput();
    normalOutput(i) = output.normalOutput;
    quadPhaseOutput_pos(i) = output.quadPhaseOutput_pos;
    quadPhaseOutput_neg(i) = output.quadPhaseOutput_neg;
    normalOutputoffset(i) = lfo.plusOffet(1/6);
end

figure;
subplot(4,1,1);
plot(normalOutput);
title("原始波形")

subplot(4,1,2);
plot(quadPhaseOutput_pos);
title("相位+90度")

subplot(4,1,3);
plot(quadPhaseOutput_neg);
title("相位-90度")

subplot(4,1,4);
plot(normalOutputoffset);
title("相位+60度")

%% SimpleDelay模块测试
[audio, fs] = audioread("data/mono_data.wav");

simpledelay = SimpleDelay();
simpledelay.reset(fs);
params = simpledelay.getParameters();
params.delayTime_mSec = 500;
simpledelay.setParameters(params);

% 测试结构体
simpledelay_output = zeros(length(audio), 1);
for i = 1 : length(audio)
    yn = simpledelay.readDelay();
    simpledelay.writeDelay(audio(i));
    simpledelay_output(i) = 0.2 * yn;
end

audiowrite('data/mono_simpledelay_output.wav', simpledelay_output, fs);

%% AudioDelay模块测试
% 2021/3/22通过测试，发现CircularBuffer中将index变为int16会导致缓存区写入出现问题
% 读取测试数据
[audio, fs] = audioread("data/delay_original.wav");

% 测试结构体标准延时
audiodelay1 = AudioDelay();
audiodelay1.reset(fs);
params = audiodelay1.getParameters();
params.feedback_Pct = 92.0;
params.leftDelay_mSec = 181;
params.rightDelay_mSec = 181.0;
params.wetLevel_dB = -6;
params.dryLevel_dB = 0;
params.interpolate = "kCubic";
audiodelay1.setParameters(params);

mono_output = zeros(length(audio), 2);
for i = 1 : length(audio)
     mono_output(i,:) = audiodelay1.processAudioFrame(audio(i,:),2,2);
end

audiowrite('data/delay_myex33.wav', mono_output, fs);

% 测试结构体ping-pong延时
% 修改参数
% audiodelay2 = AudioDelay();
% audiodelay2.reset(fs);
% params = audiodelay2.getParameters();
% params.feedback_Pct = 75.0;
% params.leftDelay_mSec = 200;
% params.rightDelay_mSec = 400.0;
% params.wetLevel_dB = -6;
% params.dryLevel_dB = 0;
% params.algorithm = 'kPingPong';
% params.reverseChannels = true;
% audiodelay2.setParameters(params);
% 
% stereo_output = zeros(length(audio), 2);
% for i = 1 : length(audio)
%     stereo_output(i,:) = audiodelay2.processAudioFrame(audio(i,:), 2, 2);
% end
% 
% audiowrite('data/pingpong_myex2.wav', stereo_output, fs);

%% ModulatedDelay模块测试
% 读取测试数据
[audio, fs] = audioread("data/chorus_original.wav");
audio = audio(:,1);
% 测试结构体单输入单输出
modulateddelay = ModulatedDelay();
modulateddelay.reset(fs);
params = modulateddelay.getParameters();
params.algorithm = "kChorus";
params.lfoRate_Hz = 0.6;
params.lfoDepth_Pct = 100;
params.feedback_Pct = 0;
params.waveform = "kSin";
modulateddelay.setParameters(params);

mono_output = zeros(length(audio), 1);
for i = 1 : length(audio)
     mono_output(i) = modulateddelay.processAudioSample(audio(i));
end

audiowrite('data/chorus_myex2.wav', mono_output, fs);

%% Flanger模块测试
[audio, fs] = audioread("data/flanger_data.wav");
%audio = audio(:,1);
% 测试结构体单输入单输出
flanger = Flanger();
flanger.reset(fs);
params = flanger.getParameters();
params.lfoRate_Hz = 0.8;
params.waveform = "kSin";
params.minDelay_mSec = 1.5;
params.modDepth_mSec = 10;
params.feedback_Pct = 50;
params.interpolate = "kCubic";
params.stereo = true;
flanger.setParameters(params);

mono_output = zeros(length(audio), 2);
for i = 1 : length(audio)
     mono_output(i,:) = flanger.processAudioFrame(audio(i,:),2,2);
end

audiowrite('data/flanger_myex23.wav', mono_output, fs);

%% Vibrato模块测试
[audio, fs] = audioread("data/vibrato_original.wav");
%audio = audio(:,1);
% 测试结构体单输入单输出
vibrato = Vibrato();
vibrato.reset(fs);
params = vibrato.getParameters();
params.lfoRate_Hz = 3.5;
params.waveform = "kTriangle";
params.modDepth_mSec = 10;
params.interpolate = "kCubic";
vibrato.setParameters(params);

mono_output = zeros(length(audio), 2);
for i = 1 : length(audio)
     mono_output(i,:) = vibrato.processAudioFrame(audio(i,:),2,2);
end

audiowrite('data/vibrato_myex33.wav', mono_output, fs);

%% Chorus模块测试
[audio, fs] = audioread("data/chorus_original.wav");
%audio = audio(:,1);
% 测试结构体单输入单输出
chorus = Chorus();
chorus.reset(fs);
params = chorus.getParameters();
params.lfoRate_Hz = 0.6;
params.waveform = "kSin";
params.minDelay_mSec = 30;
params.modDepth_mSec = 15;
params.interpolate = "kCubic";
params.stereo = false;
chorus.setParameters(params);

mono_output = zeros(length(audio), 2);
for i = 1 : length(audio)
     mono_output(i,:) = chorus.processAudioFrame(audio(i,:),2,2);
end

audiowrite('data/chorus_myex22.wav', mono_output, fs);

%% AudioDetector模块测试
% 读取测试数据
[audio, fs] = audioread("mono_data.wav");

% 测试类
audiodetector = AudioDetector();
audiodetector.reset(fs);
params = audiodetector.getParameters();
params.attackTime_mSec = 20;
params.releaseTime_mSec = 10;
params.detectMode = 0;
params.clampToUnityMax = false;
audiodetector.setParameters(params);

envelop_peak = zeros(length(audio), 1);
for i = 1 : length(audio)
     envelop_peak(i) = audiodetector.processAudioSample(audio(i));
end

audiodetector.reset(fs);
params.detectMode = 1;
audiodetector.setParameters(params);
envelop_ms = zeros(length(audio), 1);
for i = 1 : length(audio)
     envelop_ms(i) = audiodetector.processAudioSample(audio(i));
end

audiodetector.reset(fs);
params.detectMode = 2;
audiodetector.setParameters(params);
envelop_rms = zeros(length(audio), 1);
for i = 1 : length(audio)
     envelop_rms(i) = audiodetector.processAudioSample(audio(i));
end


x = 1:length(audio);
t = x / fs;

figure;
subplot(3, 1, 1);
plot(t, audio);
hold on;
plot(t, envelop_peak);

subplot(3, 1, 2);
%plot(t, audio);
hold on;
plot(t, envelop_ms);

subplot(3, 1, 3);
plot(t, audio);
hold on;
plot(t, envelop_rms);

%% CombFilter模块测试
% 读取测试数据
[audio, fs] = audioread("data/mono_data.wav");

% 测试类
combfilter = CombFilter();
combfilter.reset(fs);
params = combfilter.getParameters();
params.delayTime_mSec = 16;
params.RT60Time_mSec = 1000;
params.enableLPF = false;
params.lpf_g = 0.5;
combfilter.setParameters(params);

mono_output = zeros(length(audio), 1);
for i = 1 : length(audio)
     mono_output(i) = combfilter.processAudioSample(audio(i));
end

audiowrite('data/mono_combfilter_output.wav', mono_output, fs);

%% DelayAPF模块测试
% 读取测试数据
[audio, fs] = audioread("mono_data.wav");

% 测试类
delayapf = DelayAPF();
delayapf.reset(fs);
params = delayapf.getParameters();
params.delayTime_mSec = 16;
params.apf_g = 0.5;
params.enableLFO = true;
params.lfoRate_Hz = 100;
params.lfoDepth = 10;
params.lfoMaxModulation_mSec = 10;
delayapf.setParameters(params);

mono_output = zeros(length(audio), 1);
for i = 1 : length(audio)
     mono_output(i) = delayapf.processAudioSample(audio(i));
end

audiowrite('mono_delayapf_output.wav', mono_output, fs);

%% ImpulseConvolver模块测试
% 读取测试数据
[audio, fs] = audioread("data/mono_data.wav");

% 测试类
impulseconvolver = ImpulseConvolver();
irArray = [0.125 0.125 0.125 0.125 0.125 0.125 0.125 0.125];
impulseconvolver.setImpulseResponse(irArray, length(irArray));


mono_output = zeros(length(audio), 1);
for i = 1 : length(audio)
     mono_output(i) = impulseconvolver.processAudioSample(audio(i));
end

audiowrite('data/mono_impulse_output.wav', mono_output, fs);

%% ReverbTank模块测试
[audio, fs] = audioread("data/0.wav");

% 测试类
reverbtank = ReverbTank();
reverbtank.reset(fs);
params = reverbtank.getParameters();
params.density = "kThick"; % kThick  kThin
params.apfDelayMax_mSec = 33;
params.apfDelayWeight_Pct = 85;
params.fixeDelayMax_mSec = 81;
params.fixeDelayWeight_Pct = 100;
params.preDelayTime_mSec = 25;
params.lpf_g = 0.2;
params.kRT = 0.9;
params.lowShelf_fc = 150;
params.lowShelfBoostCut_dB = -20;
params.highShelf_fc = 4000;
params.highShelfBoostCut_dB = -6; 
params.wetLevel_dB = -12;
params.dryLevel_dB = 0;

reverbtank.setParameters(params);

mono_output = zeros(length(audio), 1);
for i = 1 : length(audio)
     mono_output(i) = reverbtank.processAudioSample(32767*audio(i));
end

audiowrite('data/0_reverb_output2.wav', mono_output/32767, fs);