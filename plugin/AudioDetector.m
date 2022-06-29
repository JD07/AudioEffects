% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% 语音包络检测模块
% 使用方法
% （1）调用构造函数初始化类的实例
% （2）调用reset方法设置采样率
% （3）调用getParameters方法，返回可调整参数结构体
% （4）设置可调整参数结构体，并利用setParameters方法设置参数
classdef AudioDetector < handle
    properties(Constant, Access = private)
        TLD_AUDIO_DETECT_MODE_PEAK = 0; % 检测峰值
        TLD_AUDIO_DETECT_MODE_MS = 1; % 检测能量
        TLD_AUDIO_DETECT_MODE_RMS = 2; % 检测均方差
    end
    
    properties(Access = public)
        audioDetectorParameters % 可调整参数结构体
        % --- attackTime_mSec % 触发时间，单位ms
        % --- releaseTime_mSec % 施放时间，单位ms
        % --- detectMode % 检测模式： 0 - peak 1 - MS 2 - RMS
        % --- detect_dB % 在dB域上检测，默认false
        % --- clampToUnityMax % 将输出钳制到1.0，当detector工作在log域上时置为false
        attackTime % 触发时间系数
        releaseTime % 施放时间系数
        sampleRate % 采样率
        lastEnvelope % 输出寄存器
    end
    
    methods(Access = public)
        % === 构造函数 ===
        function obj = AudioDetector()
            field1 = "attackTime_mSec";
            value1 = 0.0;
            field2 = "releaseTime_mSec";
            value2 = 0.0;
            field3 = "detectMode";
            value3 = 0;
            field4 = "detect_dB";
            value4 = false;
            field5 = "clampToUnityMax";
            value5 = true;
            obj.audioDetectorParameters = struct(field1,value1,field2,value2,...
                field3,value3,field4,value4,field5,value5);
            obj.attackTime = 0.0;
            obj.releaseTime = 0.0;
            obj.sampleRate = 44100;
            obj.lastEnvelope = 0;
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            obj.setSampleRate(sampleRate_);
            obj.lastEnvelope = 0.0;
        end
        
        % === 获取参数 ===
        function audioDetectorParameters = getParameters(obj)
            audioDetectorParameters = obj.audioDetectorParameters;
        end
        
        % === 参数设置 ===
        function setParameters(obj, parameters)
            obj.audioDetectorParameters = parameters;
            % 结构体参数更新
            obj.setAttackTime(obj.audioDetectorParameters.attackTime_mSec, true);
            obj.setReleaseTime(obj.audioDetectorParameters.releaseTime_mSec, true);
        end
        
        % === 采样率设置函数 ===
        function setSampleRate(obj, sampleRate_)
            if obj.sampleRate == sampleRate_
                return;
            end
            obj.sampleRate = sampleRate_;
            
            % 重新计算RC时间常数
            obj.setAttackTime(obj.audioDetectorParameters.attackTime_mSec, true);
            obj.setReleaseTime(obj.audioDetectorParameters.releaseTime_mSec, true);
        end
        
        % === 处理函数 ===
        function currEnvelope = processAudioSample(obj, xn)
            input = abs(xn);
            if obj.TLD_AUDIO_DETECT_MODE_MS == obj.audioDetectorParameters.detectMode...
                    || obj.TLD_AUDIO_DETECT_MODE_RMS == obj.audioDetectorParameters.detectMode
                input = input^2;
            end
            
            % 使用触发、释放系数进行包络跟踪
            if input > obj.lastEnvelope
                currEnvelope = obj.attackTime * obj.lastEnvelope + (1 - obj.attackTime) * input;
                %currEnvelope = obj.attackTime * (obj.lastEnvelope - input) + input;
            else
                currEnvelope = obj.releaseTime * obj.lastEnvelope + (1 - obj.releaseTime) * input;
                %currEnvelope = obj.releaseTime * (obj.lastEnvelope - input) + input;
            end
            
            if (obj.audioDetectorParameters.clampToUnityMax)
                currEnvelope = min(currEnvelope, 1.0);
            end
            
            currEnvelope = max(currEnvelope, 0.0);
            
            obj.lastEnvelope = currEnvelope;
            
            % 如果追踪的是RMS，则需要在这里开方
            if obj.TLD_AUDIO_DETECT_MODE_RMS == obj.audioDetectorParameters.detectMode
                currEnvelope = sqrt(currEnvelope);
            end
            
            % 如果不是dB，则可以返回结果了
            if ~obj.audioDetectorParameters.detect_dB
                return;
            end
            
            % 否则需要转为dB
            if currEnvelope <= 0
                currEnvelope = -96.0;
                return;
            end
            
            currEnvelope = 20 * log10(currEnvelope);
        end
        
    end
    
    methods(Access = protected)
        function setAttackTime(obj, attack_in_ms, forceCalc)
            if ~forceCalc && obj.audioDetectorParameters.attackTime_mSec == attack_in_ms
                return;
            end
            obj.audioDetectorParameters.attackTime_mSec = attack_in_ms;
            obj.attackTime = exp(-1/(attack_in_ms * obj.sampleRate * 0.001));
        end
        
        function setReleaseTime(obj, release_in_ms, forceCalc)
            if ~forceCalc && obj.audioDetectorParameters.releaseTime_mSec == release_in_ms
                return;
            end
            obj.audioDetectorParameters.releaseTime_mSec = release_in_ms;
            obj.releaseTime = exp(-1/(release_in_ms * obj.sampleRate * 0.001));
        end
    end
    
end

