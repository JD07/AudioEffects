% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% ClassATubePre类

classdef ClassATubePre < handle
    properties(Constant, Access = private)
        NUM_TUBES = 4;
    end
    
    properties(Access = public)
        parameters % 可调制参数结构体
        % --- inputLevel_dB 
        % --- saturation
        % --- asymmetry
        % --- outputLevel_dB
        % --- lowShelf_fc
        % --- lowShelfBoostCut_dB
        % --- highShelf_fc
        % --- highShelfBoostCut_dB
        triodes % TriodeClassA类数组
        shelvingFilter % TwoBandShelvingFilter类
        
        inputLevel
        outputLevel
    end

    methods(Access = public)
        % === 构造函数 ===
        function obj = ClassATubePre()
            field1 = "inputLevel_dB";
            value1 = 0.0;
            field2 = "saturation";
            value2 = 0.0;
            field3 = "asymmetry";
            value3 = 0.0;
            field4 = "outputLevel_dB";
            value4 = 0.0;
            field5 = "lowShelf_fc";
            value5 = 0.0;
            field6 = "lowShelfBoostCut_dB";
            value6 = 0.0;
            field7 = "highShelf_fc";
            value7 = 0.0;
            field8 = "highShelfBoostCut_dB";
            value8 = 0.0;
            obj.parameters = struct(field1,value1,field2,value2,field3,value3,field4,value4,...
                                    field5,value5,field6,value6,field7,value7,field8,value8);
            obj.triodes = zeros(1, obj.NUM_TUBES);
            for i = 1 : obj.NUM_TUBES
                obj.triodes(i) = TriodeClassA();
            end
            obj.shelvingFilter = TwoBandShelvingFilter();
            
            obj.inputLevel = 1.0;
            obj.outputLevel = 1.0;
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            tubeParams = obj.triodes(1).getParameters();
            tubeParams.invertOutput = true;
            tubeParams.enableHPF = true;
            tubeParams.outputGain = 1.0;
            tubeParams.saturation = 1.0;
            tubeParams.asymmetry = 0.0;
            tubeParams.enableLSF = true;
            tubeParams.lsf_Fshelf = 88.0;
            tubeParams.lsf_BoostCut_dB = -12.0;
            tubeParams.waveshaper = "kFuzzAsym";
            
            for i = 1 : obj.NUM_TUBES
               obj.triodes(i).reset(sampleRate_);
               obj.triodes(i).setParameters(tubeParams);
            end
            
            obj.shelvingFilter.reset(sampleRate_);
            
        end
        
        % === 获取参数 ===
        function parameters = getParameters(obj)
            parameters = obj.parameters;
        end
        
        % === 设置参数 ===
        function setParameters(obj, params)
            if params.inputLevel_dB ~= obj.parameters.inputLevel_dB
                obj.inputLevel = 10 ^ (params.inputLevel_dB / 20.0);
            end
            if params.outputLevel_dB ~= obj.parameters.outputLevel_dB
                obj.outputLevel = 10 ^ (params.outputLevel_dB / 20.0);
            end
            
            obj.parameters = params;
            
            sfParams = obj.shelvingFilter.getParameters();
            sfParams.lowShelf_fc = obj.parameters.lowShelf_fc;
            sfParams.lowShelfBoostCut_dB = obj.parameters.lowShelfBoostCut_dB;
            sfParams.highShelf_fc = obj.parameters.highShelf_fc;
            sfParams.highShelfBoostCut_dB = obj.parameters.highShelfBoostCut_dB;
            obj.shelvingFilter.setParameters(sfParams);
            
            tubeParams = obj.triodes(1).getParameters();
            tubeParams.saturation = obj.parameters.saturation;
            tubeParams.asymmetry = obj.parameters.asymmetry;
            
            for i = 1 : obj.NUM_TUBES
                obj.triodes(i).setParameters(tubeParams);
            end
            
        end
        
        % === 数据处理函数 ===
        function yn = processAudioSample(obj, xn)
            output1 = obj.triodes(1).processAudioSample(xn * obj.inputLevel);
            output2 = obj.triodes(2).processAudioSample(output1);
            output3 = obj.triodes(3).processAudioSample(output2);
            
            outputEQ = obj.shelvingFilter.processAudioSample(output3);
            output4 = obj.triodes(4).processAudioSample(outputEQ);
            
            yn = output4 * obj.outputLevel;
        end
        
    end
end