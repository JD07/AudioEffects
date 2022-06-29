% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% 线性缓存类
% 2021/3/10 将index调整至与C语言一致，即从0开始计算

classdef LinearBuffer < handle
    properties(Access = private)
        buffer % 数据缓存区
        bufferLength % 缓存区长度
    end
    
    methods(Access = public)
        % === 构造函数 ===
        function obj = LinearBuffer()
            obj.bufferLength = 1024;
            obj.buffer = [];
        end
        
        % === 缓存清空函数 ===
        function obj = flushBuffer(obj)
            obj.buffer(:) = 0;
        end
        
        % === 缓存创造函数 ===
        function createLinearBuffer(obj, bufferLength_)
            obj.bufferLength = bufferLength_;
            obj.buffer = zeros(1, obj.bufferLength);
            obj.flushBuffer();
        end
        
        % === 写缓存函数 ===
        % 注意，我们这里与C的习惯对齐，即index从0开始
        function writeBuffer(obj, index, input)
            if index >= obj.bufferLength
                return;
            end
            
            obj.buffer(index + 1) = input;
        end
        
        % 都缓存函数
        function output = readBuffer(obj, index)
            if index >= obj.bufferLength
                output = 0.0;
                return;
            end
            
            output = obj.buffer(index + 1);
        end
        
    end
    
end