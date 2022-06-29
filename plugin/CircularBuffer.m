% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% CircularBuffer类，实现环形缓存区域
% 特点：令缓存长度为2的n次方，这样就可以利用位操作运算实现序号范围的限制
% 注意，初始时数据缓存区为空数组，需要调用createCircularBuffer来初始化数据缓存区
% 在原代码的基础上，将整个模块改为支持cubic插值
% 目前缓存支持四种插值方法，默认为kNearest，即最近邻算法
% 插值算法中Cubic使用了x[n-1] x[n] x[n+1] x[n+2]四个值，因此，当延时小于1的时候，我们使用Square进行代替

classdef CircularBuffer < handle
    properties(Access = public)
        writeIndex  % 当前缓存位置指示器
        bufferLength % 缓存区长度
        wrapMask % 位与运算对象，确保指示器不会访问越界
        interpolate % 插值方法： kNearest, kLinear, kSquare, kCubic
        buffer % 数据缓存区
    end
    
    methods(Access = public)
        % === 构造函数 ===
        function obj = CircularBuffer()
            obj.writeIndex = 0;
            obj.bufferLength = 8;
            obj.wrapMask = obj.bufferLength - 1;
            obj.interpolate = "kNearest";
            obj.buffer = [];
        end
        
        % === 清空函数 ===
        function flushBuffer(obj)
            obj.writeIndex = 0;
            obj.buffer(:) = 0;
        end
        
        % === 重新设置缓存区长度并分配空间 ===
        % 根据输入值，从2的n次方中找出最接近的值
        function createCircularBuffer(obj, bufferLength_)
            obj.bufferLength = 2 ^ (ceil(log(bufferLength_)/log(2)));

            obj.writeIndex = 0;
            obj.wrapMask = obj.bufferLength - 1;
            obj.buffer = zeros(1, obj.bufferLength);
            obj.flushBuffer();
        end
        
        % === 向缓存中写入新的数据 ===
        function writeBuffer(obj, input)
            obj.buffer(obj.writeIndex + 1) = input;
            obj.writeIndex = obj.writeIndex + 1;
            obj.writeIndex = bitand(int32(obj.writeIndex), int32(obj.wrapMask));
        end
        
        % === 读取数据，delayInSamples为数据读取的偏移 ===
        % 若delayInSamples不为整数，则该函数根据前后位置进行插值
        function output = readBuffer(obj, delayInSamples)
            if delayInSamples < 1.0 || delayInSamples > obj.bufferLength
               error("delayInSamples is wrong")
            end
            if "kNearest" == obj.interpolate
                closestIndex = floor(delayInSamples + 0.5);
                readIndex = (obj.writeIndex) - closestIndex;
                readIndex = bitand(int32(readIndex), int32(obj.wrapMask));
                output = obj.buffer(readIndex + 1);
            elseif "kLinear" == obj.interpolate
                % 除此之外的所有值采取线性插值
                offset = floor(delayInSamples) + 1;
                fraction = offset - delayInSamples;
                readIndex = (obj.writeIndex) - offset;
                readIndex = bitand(int32(readIndex), int32(obj.wrapMask));
                previousSample = obj.buffer(readIndex + 1);
                
                readIndex = (obj.writeIndex) - (offset - 1);
                readIndex = bitand(int32(readIndex), int32(obj.wrapMask));
                nextSample = obj.buffer(readIndex + 1);

                output = circularbuffer_LinearInterpolation(previousSample, nextSample, fraction);
            elseif "kCubic" == obj.interpolate && delayInSamples >= 2
                % cubic插值需要用到x[n-1] x[n] x[n+1] x[n+2]四个值
                % 当延时小于1的时候，x[n+2]的值无法保证，所以此时不执行cubic插值
                % 计算需要读取的数据的位置
                offset = floor(delayInSamples) + 1;
                readIndex = (obj.writeIndex) - offset;
                readIndex = bitand(int32(readIndex), int32(obj.wrapMask));
                sample1 = obj.buffer(readIndex + 1);
                
                readIndex = (obj.writeIndex) - (offset - 1);
                readIndex = bitand(int32(readIndex), int32(obj.wrapMask));
                sample2 = obj.buffer(readIndex + 1);
                
                readIndex = (obj.writeIndex) - (offset - 2);
                readIndex = bitand(int32(readIndex), int32(obj.wrapMask));
                sample3 = obj.buffer(readIndex + 1);
                
                readIndex = (obj.writeIndex) - (offset + 1);
                readIndex = bitand(int32(readIndex), int32(obj.wrapMask));
                sample0 = obj.buffer(readIndex + 1);
                % fraction是插值位置距离x[n]的距离
                fraction = offset - delayInSamples;
                output = circularbuffer_CubicInterpolation(sample0, sample1, sample2, sample3, fraction);
            else
                % 除此之外的所有值采取平方插值
                offset = floor(delayInSamples) + 1;
                readIndex = (obj.writeIndex) - offset;
                readIndex = bitand(int32(readIndex), int32(obj.wrapMask));
                sample1 = obj.buffer(readIndex + 1);
                
                readIndex = (obj.writeIndex) - (offset - 1);
                readIndex = bitand(int32(readIndex), int32(obj.wrapMask));
                sample2 = obj.buffer(readIndex + 1);
                
                readIndex = (obj.writeIndex) - (offset + 1);
                readIndex = bitand(int32(readIndex), int32(obj.wrapMask));
                sample0 = obj.buffer(readIndex + 1);
                
                % fraction是插值位置距离x[n]的距离
                fraction = offset - delayInSamples;

                output = circularbuffer_SquareInterpolation(sample0, sample1, sample2, fraction);
            end

        end
        
        % === 插值功能控制 ===
        function setInterpolate(obj, b)
            obj.interpolate = b;
        end
        
    end
end

function y = circularbuffer_LinearInterpolation(x0, x1, fraction)
    if fraction < 0.0 || fraction > 1.0
        error("ErrorXXX: fraction should in range[0 1]!\n");
    end
    y = fraction * x1 + (1 - fraction) * x0;
    end

function y = circularbuffer_SquareInterpolation(x_1, x0, x1, fraction)
    if fraction < 0.0 || fraction > 1.0
        error("ErrorXXX: fraction should in range[0 1]!\n");
    end

    frsq = fraction * fraction;
    a0 = x0;
    a1 = 0.5 * (x1 - x_1);
    a2 = 0.5 * (x1 - 2*x0 + x_1);

    y = a2 * frsq + a1 * fraction + a0;
end

function y = circularbuffer_CubicInterpolation(x_1, x0, x1, x2, fraction)
    if fraction < 0.0 || fraction > 1.0
        error("ErrorXXX: fraction should in range[0 1]!\n");
    end

    frsq = fraction * fraction;
    a0 = -0.5 * x_1 + 1.5 * x0 - 1.5 * x1 + 0.5 * x2;
    a1 = x_1 - 2.5 * x0 + 2 * x1 - 0.5 * x2;
    a2 = -0.5 * x_1 + 0.5 * x1;
    a3 = x0;

    y = a0 * fraction * frsq + a1 * frsq + a2 * fraction + a3; 
end
