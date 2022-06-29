% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% DynamicsProcessor

classdef DynamicsProcessor < handle
    properties(Access = public)
        parameters % 可调制参数结构体
        % --- ratio
        % --- threshold_dB
        % --- kneeWidth_dB
        % --- hardLimitGate
        % --- softKnee
        % --- enableSidechain
        % --- calculation
        % --- attackTime_mSec
        % --- releaseTime_mSec
        % --- outputGain_dB
        % --- gainReduction
        % --- gainReduction_dB
        detector % AudioDetector类
        sidechainInputSample
    end
    
    
    % === 构造函数 ===
    methods(Access = public)
        % === 构造函数 ===
        function obj = DynamicsProcessor()
            % 参数结构体默认参数
            obj.parameters = struct();
            obj.parameters.ratio = 50.0;
            obj.parameters.threshold_dB = -10.0;
            obj.parameters.kneeWidth_dB = 10.0;
            obj.parameters.hardLimitGate = false;
            obj.parameters.softKnee = true;
            obj.parameters.enableSidechain = false;
            obj.parameters.calculation = "kCompressor";
            obj.parameters.attackTime_mSec = 0.0;
            obj.parameters.releaseTime_mSec = 0.0;
            obj.parameters.outputGain_dB = 0.0;
            obj.parameters.gainReduction = 1.0;
            obj.parameters.gainReduction_dB = 0.0;
            
            obj.detector = AudioDetector();
            
            obj.sidechainInputSample = 0.0;
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            obj.sidechainInputSample = 0.0;
            obj.detector.reset(sampleRate_);
            detectorParams = obj.detector.getParameters();
            detectorParams.clampToUnityMax = false;
            detectorParams.detect_dB = true;
            obj.detector.setParameters(detectorParams);
        end
       
        % === 使能旁路输入 ===
        function enableAuxInput(obj, enableAuxInput)
            obj.parameters.enableSidechain = enableAuxInput;
        end
        
        % === 处理旁路输入 ===
        function yn = processAuxInputAudioSample(obj, xn)
            obj.sidechainInputSample = xn;
            yn = xn;
        end
        
        % === 获取参数 ===
        function parameters = getParameters(obj)
            parameters = obj.parameters;
        end
        
        % === 设置参数 ===
        function setParameters(obj, parameters_)
            obj.parameters = parameters_;
            
            detectorParams = obj.detector.getParameters();
            detectorParams.attackTime_mSec = obj.parameters.attackTime_mSec;
            detectorParams.releaseTime_mSec = obj.parameters.releaseTime_mSec;
            obj.detector.setParameters(detectorParams);
        end
        
        % === 数据处理函数 ===
        function yn = processAudioSample(obj, xn)
            if obj.parameters.enableSidechain
                detect_dB = obj.detector.processAudioSample(obj.sidechainInputSample);
            else
                detect_dB = obj.detector.processAudioSample(xn);
            end
            
            gr = obj.computeGain(detect_dB);
            
            makeupGain = 10.0^(obj.parameters.outputGain_dB / 20.0);
            
            yn = xn * gr * makeupGain;
        end
        
        % === 增益计算函数 ===
        function gain = computeGain(obj, detect_dB)
            parameters_ = obj.parameters;
            if "kCompressor" == parameters_.calculation
                if ~parameters_.softKnee
                    if detect_dB <= parameters_.threshold_dB
                        output_dB = detect_dB; 
                    else
                        if parameters_.hardLimitGate
                            output_dB = parameters_.threshold_dB;
                        else
                            output_dB = parameters_.threshold_dB +... 
                                        (detect_dB - parameters_.threshold_dB) / parameters_.ratio;
                        end
                    end
                else
                    if 2.0 * (detect_dB - parameters_.threshold_dB) < -parameters_.kneeWidth_dB
                        output_dB = detect_dB;
                    elseif 2.0 * abs(detect_dB - parameters_.threshold_dB) <= parameters_.kneeWidth_dB
                        if parameters_.hardLimitGate
                            output_dB = detect_dB -... 
                                (detect_dB - parameters_.threshold_dB + (parameters_.kneeWidth_dB / 2.0))^2 / (2.0*parameters_.kneeWidth_dB);
                        else
                            output_dB = detect_dB +...
                                (((1.0 / parameters_.ratio) - 1.0) * (detect_dB - parameters_.threshold_dB + (parameters_.kneeWidth_dB / 2.0))^2.0) / (2.0*parameters_.kneeWidth_dB);
                        end
                    elseif 2.0 * (detect_dB - parameters_.threshold_dB) > parameters_.kneeWidth_dB
                        if parameters_.hardLimitGate
                            output_dB = parameters_.threshold_dB;
                        else
                            output_dB = parameters_.threshold_dB +...
                                        (detect_dB - parameters_.threshold_dB) /  parameters_.ratio;
                        end
                    end
                end
            elseif "kDownwardExpander" == parameters_.calculation
                if ~parameters_.softKnee || parameters_.hardLimitGate
                    if detect_dB >= parameters_.threshold_dB
                        output_dB = detect_dB;
                    else
                        if parameters_.hardLimitGate
                            output_dB = -1.0e34;
                        else
                            output_dB = parameters_.threshold_dB +...
                                (detect_dB - parameters_.threshold_dB)*parameters_.ratio;
                        end
                    end
                else
                    if 2.0 * (detetc_dB - parameters_.threshold_dB) > parameters_.kneeWidth_dB
                        output_dB = detect_dB;
                    elseif 2.0 * abs(detect_dB - parameters_.threshold_dB) > -parameters_.kneeWidth_dB
                        output_dB = (parameters_.ratio - 1) *...
                                (detect_dB - parameters_.threshold_dB - (parameters_.kneeWidth_dB / 2.0))^2.0 / (2.0*parameters_.kneeWidth_dB);
                    elseif 2.0 * (detect-dB - parameters_.threshold_dB) <= -parameters_.kneeWidth_dB
                        output_dB = parameters_.threshold_dB + (detect_dB - parameters_.threshold_dB)*parameters_.ratio;
                    end
                end
            end
            
            parameters_.gainReduction_dB = output_dB - detect_dB;
            parameters_.gainReduction = 10 ^ (parameters_.gainReduction_dB / 20.0);
            
            obj.parameters = parameters_;
            
            gain = obj.parameters.gainReduction;
        end
        
    end
end


