% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% AnalogFIRFilter模块

classdef AnalogFIRFilter < handle
    properties(Constant, Access = private)
        IR_LEN = 512;
    end
    
    properties(Access = public)
        parameters % 可调制参数结构体
        % --- filterType % kLPF1, kHPF1, kLPF2, kHPF2, kBPF2, kBSF2
        % --- fc
        % --- Q
        convolver % ImpulseConvolver模块
        analogMagArray % 模拟幅度相应数组
        irArray % 冲激响应系数数组
        sampleRate % 采样率
    end
    
    methods
        % === 构造函数 ===
        function obj = AnalogFIRFilter()
            field1 = "filterType";
            value1 = "kLPF1";
            field2 = "fc";
            value2 = 0.0;
            field3 = "Q";
            value3 = 0.0;
            obj.parameters = struct(field1,value1,field2,value2,field3,value3);
            
            obj.convolver = ImpulseConvolver();
            obj.analogMagArray = zeros(1, obj.IR_LEN);
            obj.irArray = zeros(1, obj.IR_LEN);
            obj.sampleRate = 0.0;
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            obj.sampleRate = sampleRate_;
            obj.convolver.reset(sampleRate_);
            obj.convolver.init(obj.IR_LEN);
            
            obj.analogMagArray(:) = 0;
            obj.irArray(:) = 0;
        end
        
        % === 数据处理函数 ===
        function yn = processAudioSample(obj, xn)
            yn = obj.convolver.processAudioSample(xn);
        end
        
        % === 获取参数 ===
        function parameters = getParameters(obj)
            parameters = obj.parameters;
        end
        
        % === 参数设置 ===
        function setParameters(obj, parameters_)
            if parameters_.fc ~= obj.parameters.fc ||...
               parameters_.Q ~= obj.parameters.Q ||...
               parameters_.filterType ~= obj.parameters.filterType
                
                analogFilterData.sampleRate = obj.sampleRate;
                analogFilterData.magArray = obj.analogMagArray;
                analogFilterData.dftArrayLen = obj.IR_LEN;
                analogFilterData.mirrorMag = false;
                
                analogFilterData.filterType = parameters_.filterType;
                analogFilterData.fc = parameters_.fc;
                analogFilterData.Q = parameters_.Q;
                
                calculateAnalogMagArray(analogFilterData);
                freqSample(obj.IR_LEN, obj.analogMagArray, obj.irArray, 1);
                
                obj.convolver.setImpulseResponse(obj.irArray, obj.IR_LEN);
            end
            
            obj.parameters = parameters_;
        end
        
    end
    
end

function magData = calculateAnalogMagArray(magData)
end

function freqSample(N, A, h, symm)
end
