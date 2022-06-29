% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% 调制音频延时模块
% 实现三种音效效果：flanger,chorus,vibrato
% 使用方法：
% （1）调用构造函数初始化类的实例
% （2）调用reset方法设置采样率
% （3）调用getParameters方法，返回可调整参数结构体
% （4）设置可调整参数结构体，并利用setParameters方法设置参数

% 2021/2/23 通过单声道音频测试

classdef ModulatedDelay < handle
    properties(Access = public)
        parameters % 可调整参数结构体
        % --- algorithm % kFlanger, kChorus, kVibrato
        % --- lfoRate_Hz
        % --- lfoDepth_Pct
        % --- feedback_Pct
        % --- waveform % LFO波形 kSin, kTriangle, kSaw
        delay % AudioDelay类
        lfo % LFO类
    end
    
    methods(Access = public)
        % === 构造函数 ===
        function obj = ModulatedDelay()
            obj.parameters = struct();
            obj.parameters.algorithm = "kFlanger";
            obj.parameters.lfoRate_Hz = 0.0;
            obj.parameters.lfoDepth_Pct = 0.0;
            obj.parameters.feedback_Pct = 0.0;
            obj.parameters.waveform = "kTriangle";
            
            obj.delay = AudioDelay();
            obj.lfo = LFO();
        end
        
        % === 重置函数 ===
        function obj = reset(obj, sampleRate_)
            % 创建新的延迟缓存区，长度为100ms
            obj.delay.reset(sampleRate_);
            % flanger、chorus以及vibrato的延时一般都不会超过30ms（超过就会像echo）
            obj.delay.createDelayBuffers(sampleRate_, 100);
            
            % lfo
            obj.lfo.reset(sampleRate_);
            params = obj.lfo.getParameters();
            params.waveform = "kTriangle";
            obj.lfo.setParameters(params);
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
            adParams.feedback_Pct = obj.parameters.feedback_Pct;
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
            minDelay_mSec = 0.0;
            maxDepth_mSec = 0.0;
            
            % 设置延时时间，干湿系数以及反馈系数
            if "kFlanger" == obj.parameters.algorithm
                minDelay_mSec = 2.5;
                maxDepth_mSec = 12.5;
                params.wetLevel_dB = -3.0;
                params.dryLevel_dB = -3.0;
            elseif "kChorus" == obj.parameters.algorithm
                minDelay_mSec = 30.0;
                maxDepth_mSec = 42.0;
                params.wetLevel_dB = -3.0;
                params.dryLevel_dB = -0.0;
                params.feedback_Pct = 0.0;
            elseif "kVibrato" == obj.parameters.algorithm
                minDelay_mSec = 0.0;
                maxDepth_mSec = 10;
                params.wetLevel_dB = 0.0;
                params.dryLevel_dB = -96.0;
                params.feedback_Pct = 0.0;
            end
            
            % 计算调制延时时间
            depth = obj.parameters.lfoDepth_Pct / 100;
            modulationMin = minDelay_mSec;
            modulationMax = minDelay_mSec + maxDepth_mSec;
            
            % --flanger -unipolar
            if "kFlanger" == obj.parameters.algorithm
                params.leftDelay_mSec = doUnipolarModulationFromMin(...
                    depth * bipolarToUnipolar(lfoOutput.normalOutput), modulationMin, modulationMax);
            else
                params.leftDelay_mSec = doBipolarModulation(...
                    depth * lfoOutput.normalOutput, modulationMin, modulationMax);
            end
            
            % 右声道延迟匹配左声道
            params.rightDelay_mSec = params.leftDelay_mSec;
            % 调制延时模块
            obj.delay.setParameters(params);
            % 调用延时模块的处理函数
            outputFrame = obj.delay.processAudioFrame(inputFrame, inputChannels, outputChannels);
        end
        
    end
    
end

function output = bipolarToUnipolar(input)
    output = 0.5 * input + 0.5;
end

function output = doUnipolarModulationFromMin(unipolarModulatorValue, minValue, maxValue)
    unipolarModulatorValue = min(unipolarModulatorValue, 1.0);
    unipolarModulatorValue = max(unipolarModulatorValue, 0.0);
    
    output = unipolarModulatorValue*(maxValue - minValue) + minValue;
end

function output = doBipolarModulation(bipolarModulatorValue, minValue, maxValue)
    bipolarModulatorValue = min(bipolarModulatorValue, 1.0);
    bipolarModulatorValue = max(bipolarModulatorValue, -1.0);
    
    halfRange = (maxValue - minValue) / 2.0;
    midpoint = halfRange + minValue;

    output = bipolarModulatorValue * halfRange + midpoint;
end


