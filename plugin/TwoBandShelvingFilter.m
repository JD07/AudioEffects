% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% TwoBandShelvingFilter模块

classdef TwoBandShelvingFilter < handle
    properties(Access = public)
        parameters % 可调制参数结构体
        % --- lowShelf_fc
        % --- lowShelfBoostCut_dB
        % --- highShelf_fc
        % --- highShelfBoostCut_dB
        lowShelfFilter % AudioFilter类
        highShelfFilter % AuidoFilter类
    end
    
    methods(Access = public)
        % === 构造函数 ===
        function obj = TwoBandShelvingFilter()
            obj.parameters = struct();
            obj.parameters.lowShelf_fc = 0.0;
            obj.parameters.lowShelfBoostCut_dB = 0.0;
            obj.parameters.highShelf_fc = 0.0;
            obj.parameters.highShelfBoostCut_dB = 0.0;
            
            obj.lowShelfFilter = AudioFilter();
            obj.highShelfFilter = AudioFilter();
            
            params = obj.lowShelfFilter.getParameters();
            params.algorithm = "kLowShelf";
            obj.lowShelfFilter.setParameters(params);
            
            params = obj.highShelfFilter.getParameters();
            params.alogorithm = "kHiShelf";
            obj.highShelfFilter.setParameters(params);
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            obj.lowShelfFilter.reset(sampleRate_);
            obj.highShelfFilter.reset(sampleRate_);
        end
        
        % === 数据处理函数 ===
        function yn = processAudioSample(obj, xn)
            filteredSignal = obj.lowShelfFilter.processAudioSample(xn);
            filteredSignal = obj.highShelfFilter.processAudioSample(filteredSignal);
            
            yn = filteredSignal;
        end
        
        % === 参数获取函数 ===
        function parameters = getParameters(obj)
            parameters = obj.parameters;
        end
        
        % === 参数设置函数 ===
        function setParameters(obj, params)
            obj.parameters = params;
            filterParams = obj.lowShelfFilter.getParameters();
            filterParams.fc = obj.parameters.lowShelf_fc;
            filterParams.boostCut_dB = obj.parameters.lowShelfBoostCut_dB;
            obj.lowShelfFilter.setParameters(filterParams);
            
            filterParams = obj.highShelfFilter.getParameters();
            filterParams.fc = obj.parameters.highShelf_fc;
            filterParams.boostCut_dB = obj.parameters.highShelfBoostCut_dB;
            obj.highShelfFilter.setParameters(filterParams);
        end
        
    end
    
end