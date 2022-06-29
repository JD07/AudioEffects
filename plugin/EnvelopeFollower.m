% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% EnvelopeFollower

classdef EnvelopeFollower < handle
    properties(Constant,Access = private)
        kMaxFilterFrequency = 20480.0; % 20Hz的十倍频
    end
    
    properties(Access = public)
        parameters % 可调参数结构体
        % --- fc 
        % --- Q
        % --- attackTime_mSec
        % --- releaseTime_mSec
        % --- threshold_dB
        % --- sensitivity
        filter % ZVAFilter类
        detector % AudioDetector类
    end
    
    methods(Access = public)
        % === 构造函数 ===
        function obj = EnvelopeFollower()
            obj.parameters = struct();
            obj.parameters.fc = 0.0;
            obj.parameters.Q = 0.707;
            obj.parameters.attackTime_mSec = 10.0;
            obj.parameters.releaseTime_mSec = 10.0;
            obj.parameters.threshold_dB = 0.0;
            obj.parameters.sensitivity = 1.0;

            obj.filter = ZVAFilter();
            obj.detector = AudioDetector;
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            obj.filter.reset(sampleRate_);
            obj.detector.reset(sampleRate_);
        end
        
        % === 获取参数 ===
        function parameters = getParameters(obj)
            parameters = obj.parameters;
        end
        
        % === 设置参数 ===
        function setParameters(obj, params)
            filterParams = obj.filter.getParameters();
            adParams = obj.detector.getParameters();
            
            if params.fc ~= obj.parameters.fc || params.Q ~=obj.parameters.Q
                filterParams.fc = params.fc;
                filterParams.Q = params.Q;
                obj.filter.setParameters(filterParams);
            end
            
            if params.attackTime_mSec ~= obj.parameters.attackTime_mSec ||...
                    params.releaseTime_mSec ~= obj.parameters.releaseTime_mSec
                adParams.attackTime_mSec = params.attackTime_mSec;
                adParams.releaseTime_mSec = params.releaseTime_mSec;
                obj.detector.setParameters(adParams);
            end
            obj.parameters = params;
        end
        
        % === 数据处理 ===
        function yn = processAudioSample(obj, xn)
            threshValue = 10^(obj.parameters.threshold_dB / 20.0);
            
            detect_dB = obj.detector.processAudioSample(xn);
            detectValue = 10^(detect_dB / 20.0);
            deltaValue = detectValue - threshValue;
            
            filterParams = obj.filter.getParameters();
            filterParams.fc = obj.parameters.fc;
            
            if deltaValue > 0
                modulatorValue = deltaValue * obj.parameters.sensitivity;
                filterParams.fc = doUnipolarModulationFromMin(modulatorValue,...
                                        obj.parameters.fc, obj.kMaxFilterFrequency);
            end
            
            obj.filter.setParameters(filterParams);
            
            yn = obj.filter.processAudioSample(xn);
        end
       
        
    end
end