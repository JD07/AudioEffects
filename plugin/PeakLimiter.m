% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% PeakLimiter类

classdef PeakLimiter < handle
    properties(Constant, Access = private)
        softknee = true;
        kneeWidth_dB = 10.0;
    end
    
    properties(Access = public)
        detector % AudioDetector类
        threshold_dB
        makeUpGain_dB
    end
    
    methods(Access = public)
        % === 构造函数 ===
        function obj = PeakLimiter()
            obj.detector = AudioDetector();
            obj.threshold_dB = 0.0;
            obj.makeUpGain_dB = 0.0;
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            obj.detector.setSampleRate(sampleRate_);
            
            detectorParams = obj.detector.getParameters();
            detectorParams.detect_dB = true;
            detectorParams.attackTime_mSec = 5.0;
            detectorParams.releaseTime_mSec = 25.0;
            detectorParams.clampToUnityMax = false;
            detectorParams.detectMode = "ENVELOPE_DETECT_MODE_PEAK";
            obj.detector.setParameters(detectorParams);
        end
        
        % === 数据处理函数 ===
        function yn = processAudioSample(obj, xn)
            yn = dB2Raw(obj.makeUpGain_dB) * xn * obj.computeGain(obj.detector.processAudioSample(xn));
        end
        
        % === 增益计算函数 ===
        function gain = computeGain(obj, detect_dB)
            if ~obj.softknee
                if detect_dB <= obj.threshold_dB
                    output_dB = detect_dB;
                else
                    output_dB = obj.threshold_dB;
                end
            else
                if 2.0 * (detect_dB - obj.threshold_dB) < -obj.kneeWidth_dB
                    output_dB = detect_dB;
                elseif 2.0 * abs(detect_dB - obj.threshold_dB) <= obj.kneeWidth_dB
                    output_dB = detect_dB -...
                        (detect_dB - obj.threshold_dB + obj.kneeWidth_dB / 2.0)^2 / (2.0 * obj.kneeWidth_dB);
                elseif 2.0 * (detect_dB - obj.threshold_dB) > obj.kneeWidth_dB
                    output_dB = obj.threshold_dB;
                end
            end
            
            delta_dB = output_dB - detect_dB;
            gain = 10^(delta_dB / 20.0);
        end
        
        % === 阈值设置 ===
        function setThreshold_dB(obj, threshold_dB_)
            obj.threshold_dB = threshold_dB_;
        end
        
        % === 增益设置 ===
        function setMakeUpGain_dB(obj, makeUpGain_dB_)
            obj.makeUpGain_dB = makeUpGain_dB_;
        end
        
    end
end

function gain = dB2Raw(dB)
    gain = 10 ^ (dB / 20.0);
end
