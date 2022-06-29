% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% ImpulseConvolver类
% 该类的作用就是将冲激响应系数与输入数据进行卷积
% 2021/3/10 通过测试

classdef ImpulseConvolver < handle
    properties(Access = public)
        signalBuffer % CircularBuffer类
        irBuffer % LinearBuffer类
        length; % 长度
    end
    
    methods(Access = public)
        % === 构造函数 ===
        function obj = ImpulseConvolver()
            obj.signalBuffer = CircularBuffer();
            obj.irBuffer = LinearBuffer();
            obj.init(512);
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
           obj.signalBuffer.flushBuffer(); 
        end
        
        % === 数据处理函数 ===
        function yn = processAudioSample(obj, xn)
            obj.signalBuffer.writeBuffer(xn);
            
            output = 0;
            for i = 0 : obj.length-1
                signal = obj.signalBuffer.readBuffer(i);
                irrrrr = obj.irBuffer.readBuffer(i);
                output = output + signal * irrrrr;
            end
            
            yn = output;
        end
        
        % === 初始化函数 ===
        function init(obj, lengthPowerOfTwo)
            obj.length = lengthPowerOfTwo;
            
            obj.signalBuffer.createCircularBuffer(lengthPowerOfTwo);
            obj.irBuffer.createLinearBuffer(lengthPowerOfTwo);
        end
        
        % === 设置冲激响应 ===
        function setImpulseResponse(obj, irArray, lengthPowerOfTwo)
            if lengthPowerOfTwo ~= obj.length
                obj.length = lengthPowerOfTwo;
                obj.signalBuffer.createCircularBuffer(lengthPowerOfTwo);
                obj.irBuffer.createLinearBuffer(lengthPowerOfTwo);
            end
            
            % 装填IR缓存，注意索引进行调整以符合C语言习惯
            for i = 0 : obj.length-1
                obj.irBuffer.writeBuffer(i, irArray(i+1));
            end
        end
        
    end
end