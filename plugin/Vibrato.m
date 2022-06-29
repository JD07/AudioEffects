% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% Vibrato模块
% 使用方法：
% （1）调用构造函数初始化类的实例
% （2）调用reset方法设置采样率
% （3）调用getParameters方法，返回可调整参数结构体
% （4）设置可调整参数结构体，并利用setParameters方法设置参数

classdef Vibrato < handle
    properties(Access = public)
        parameters % 可调参数结构体
        % --- lfoRate_Hz
        % --- waveform % LFO波形 kSin, kTriangle, kSaw
        % --- modDepth_mSec; % 延时模块变动范围
        % --- interpolate % 插值方法： kNearest, kLinear, kSquare, kCubic
        delay % AudioDelay类
        lfo % LFO类
    end
    
    methods(Access = public)
        % === 构造函数 ===
        function obj = Vibrato()
            obj.parameters = struct();
            obj.parameters.lfoRate_Hz = 0.0;
            obj.parameters.waveform = "kSin";
            obj.parameters.modDepth_mSec = 0.0;
            obj.parameters.interpolate = "kNearest";
            
            obj.delay = AudioDelay();
            obj.lfo = LFO();
        end
        
        % === 重置函数 ===
        function obj = reset(obj, sampleRate_)
            obj.delay.reset(sampleRate_);
            % flanger、chorus以及vibrato的延时一般都不会超过100ms（超过就会像echo）
            obj.delay.createDelayBuffers(sampleRate_, 100);

            % lfo
            obj.lfo.reset(sampleRate_);            
        end
        
        % === 获取参数 ===
        function parameters = getParameters(obj)
            parameters = obj.parameters;
        end
        
        % === 参数设置函数 ===
        function setParameters(obj, parameters_)
            obj.parameters = parameters_;
            
            lfoParams = obj.lfo.getParameters();
            lfoParams.frequency_Hz = obj.parameters.lfoRate_Hz;
            lfoParams.waveform = obj.parameters.waveform;
            obj.lfo.setParameters(lfoParams);
            
            adParams = obj.delay.getParameters();
            adParams.feedback_Pct = 0.0;
            adParams.interpolate = obj.parameters.interpolate;
            obj.delay.setParameters(adParams);
        end
        
        % === 处理函数 ===
        function output = processAudioSample(obj, input)
            output = obj.processAudioFrame(input, 1, 1);
        end
        
        % === 音频帧处理函数 ===
        function outputFrame = processAudioFrame(obj, inputFrame, inputChannels, outputChannels)
            if inputChannels == 0 || outputChannels == 0
                error('Error : InputChannels and outputChannels should not be 0');
            end
            
            % 启动低频振荡器
            lfoOutput = obj.lfo.renderAudioOutput();
            
            % 设置延迟调制
            params = obj.delay.getParameters();
            modDepth_mSec = obj.parameters.modDepth_mSec;
            params.wetLevel_dB = -0.0;
            params.dryLevel_dB = -96.0;
            
            % 计算调制延时时间            
            params.leftDelay_mSec = 0.5 * (lfoOutput.normalOutput + 1) * modDepth_mSec;
            params.rightDelay_mSec = params.leftDelay_mSec;
            
            % 调制延时模块
            obj.delay.setParameters(params);
            % 调用延时模块的处理函数
            outputFrame = obj.delay.processAudioFrame(inputFrame, inputChannels, outputChannels);
        end
    end
end