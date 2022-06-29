% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% 标准一阶低通滤波器
% 使用方法：
% （1）调用构造函数初始化类的实例
% （2）调用getParameters方法，返回可调整参数结构体
% （3）设置可调整参数结构体，并利用setParameters方法设置参数
% （4）调用reset方法重置类的状态，完成初始化
% 2021/2/22 通过测试
classdef SimpleLPF < handle
    properties(Access = private)
        simpleLPFParameters % 可调整参数结构体
        % --- g % 系数
        state
    end
    
    methods(Access = public)
        % === 构造函数 ===
        function obj = SimpleLPF()
            obj.simpleLPFParameters.g = 0.0;
            obj.state = 0.0;
        end

        % === 重置函数 ===
        function reset(obj,~)
            obj.state = 0;
        end
        
        % === 获取参数 ===
        function params = getParameters(obj)
            params = obj.simpleLPFParameters;
        end
        
        % === 参数设置 ===
        function setParameters(obj, params)
            obj.simpleLPFParameters = params;
        end
        
        % === 模块处理数据函数 ===
        function output = processAudioSample(obj, xn)
            output = (1 - obj.simpleLPFParameters.g) * xn + obj.simpleLPFParameters.g * obj.state;
            obj.state = output;
        end
        
    end
end