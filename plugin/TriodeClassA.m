% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% TriodeClassAParameters类

classdef TriodeClassA < handle
    properties
        parameters % 可调整参数结构体
        % --- waveshaper % kSoftClip, kArcTan, kFuzzAsym
        % --- saturation
        % --- asymmetry
        % --- outputGain
        % --- invertOutput
        % --- enableHPF
        % --- enableLSF
        % --- hpf_Fc
        % --- lsf_Fshelf
        % --- lsf_BoostCut_dB
        outputHPF % AudioFilter类
        outputLSF % AudioFilter类
    end
    
    methods
        % === 构造函数 ===
        function obj = TriodeClassA()
            field1 = "waveshaper";
            value1 = "kSoftClip";
            field2 = "saturation";
            value2 = 1.0;
            field3 = "asymmetry";
            value3 = 0.0;
            field4 = "outputGain";
            value4 = 1.0;
            field5 = "invertOutput";
            value5 = true;
            field6 = "enableHPF";
            value6 = true;
            field7 = "enableLSF";
            value7 = false;
            field8 = "hpf_Fc";
            value8 = 1.0;
            field9 = "lsf_Fshelf";
            value9 = 80.0;
            field10 = "lsf_BoostCut_dB";
            value10 = 0.0;
            obj.parameters = struct(field1,value1,field2,value2,field3,value3,field4,value4,...
               field5,value5,field6,value6,field7,value7,field8,value8,field9,value9,field10,value10);
            
            obj.outputHPF = AudioFilter();
            obj.outputLSF = AudioFilter();
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            obj.outputHPF.reset(sampleRate_);
            obj.outputLSF.reset(sampleRate_);
        end
        
        % === 参数获取 ===
        function parameters = getParameters(obj)
            parameters = obj.parameters;
        end
        
        % === 参数设置 ===
        function setParameters(obj, params)
            obj.parameters = params;
            
            filterParams = obj.outputHPF.getParameters();
            filterParams.algorithm = "kHPF1";
            filterParams.fc = obj.parameters.hpf_Fc;
            obj.outputHPF.setParameters(filterParams);
            
            filterParams.algorithm = "kLowShelf";
            filterParams.fc = obj.parameters.lsf_Fshelf;
            filterParams.boostCut_dB = obj.parameters.lsf_BoostCut_dB;
            obj.outputLSF.setParameters(filterParams);
        end
        
        % === 参数处理函数 ===
        function yn = processAudioSample(obj, xn)
            if "kSoftClip" == obj.parameters.waveshaper
                output = softClipWaveShaper(xn, obj.parameters.saturation);
            elseif "kArcTan" == obj.parameters.waveshaper
                output = atanWaveShaper(xn, obj.parameters.saturation);
            elseif "kFuzzAsym" == obj.parameters.waveshaper
                output = fuzzExp1WaveShaper(xn, obj.parameters.saturation, obj.parameters.asymmetry);
            end
            
            if obj.parameters.invertOutput
                output = -output;
            end
            
            if obj.parameters.enableHPF
                output = obj.outputHPF.processAudioSample(output);
            end
            
            if obj.parameters.enableLSF
                output = obj.outputLSF.processAudioSample(output);
            end
            
            output = output * obj.parameters.outputGain;
            
            yn = output;
        end
        
    end
end