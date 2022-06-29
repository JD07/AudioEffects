% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% LRFilterBank

classdef LRFilterBank < handle
    properties(Access = public)
        parameters % 可调整参数结构体
        % --- splitFrequency
        lpFilter % AudioFilter类
        hpFilter % AudioFilter类
    end
    
    methods(Access = public)
        % === 构造函数 ===
        function obj = LRFilterBank()
            field1 = "splitFrequency";
            value1 = 1000.0;
            obj.parameters = struct(field1,value1);
            obj.lpFilter = AudioFilter();
            obj.hpFilter = AudioFilter();
            
            params = obj.lpFilter.getParameters();
            params.algorithm = "kLWRLPF2";
            obj.lpFilter.setParameters(params);
            
            params = obj.hpFilter.getParameters();
            params.algorithm = "kLWRHPF2";
            obj.hpFilter.setParameters(params);
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            obj.lpFilter.reset(sampleRate_);
            obj.hpFilter.reset(sampleRate_);
        end
        
        % === 数据处理函数 ===
        % 不做任何处理
        function yn = processAudioSample(xn)
            yn = xn;
        end
        
        % === 滤波器组处理函数 ===
        % 输出结构体，由低频和高频两部分组成
        function output = processFilterBank(obj, xn)
            output.LFOut = obj.lpFilter.processAudioSample(xn);
            % HPF的输出取反，这样后面重新组合时能保证正确的相位和幅度响应
            output.HFOut = -obj.hpFilter.processAudioSample(xn);
        end
        
        % === 获取参数 ===
        function parameters = getParameters(obj)
            parameters = obj.parameters;
        end
        
        % === 设置参数 ===
        function setParameters(obj, parameters_)
            obj.parameters = parameters_;
            
            params = obj.lpFilter.getParameters();
            params.fc = obj.parameters.splitFrequency;
            obj.lpFilter.setParameters(params);
            
            params = obj.hpFilter.getParameters();
            params.fc = obj.parameters.splitFrequency();
            obj.hpFilter.setParameters(params);
        end
        
    end
end