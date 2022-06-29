% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% 标准延时模块
% 设计思路：需要利用CircularBuffer实现缓存区，缓存长度和延迟长度根据时间（单位：ms）计算
% 使用方法
% （1）调用构造函数初始化类的实例
% （2）调用reset方法设置采样率
% （3）调用getParameters方法，返回可调整参数结构体
% （4）设置可调整参数结构体，并利用setParameters方法设置参数

% 2021/2/22通过测试

classdef SimpleDelay < handle
    properties(Access = public)
        simpleDelayParameters % 可调制参数结构体
        % --- delayTime_mSec % 延时时间长度（单位：ms）
        % --- interpolate % 插值类型 kNearest, kLinear, kSquare, kCubic
        % --- delay_Samples % 延时对应的数据点数
        sampleRate % 采样率
        samplesPerMSec % 每ms对应的数据点个数
        bufferLength_mSec % 缓存区时间长度（单位：ms）
        bufferLength % 缓存区长度（数据点个数）
        delayBuffer % 缓存区,类变量
    end
    
    methods(Access = public)
        % === 构造函数 ===
        function obj = SimpleDelay()
            obj.simpleDelayParameters = struct();
            obj.simpleDelayParameters.delayTime_mSec = 0.0;
            obj.simpleDelayParameters.interpolate = "kNearest";
            obj.simpleDelayParameters.delay_Samples = 0.0;
            
            obj.sampleRate = 0.0;
            obj.samplesPerMSec = 0.0;
            obj.bufferLength_mSec = 0.0;
            obj.bufferLength = 0;
            obj.delayBuffer = CircularBuffer();
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            if obj.sampleRate == sampleRate_
                obj.delayBuffer.flushBuffer();
            else
                obj.createDelayBuffer(sampleRate_, obj.bufferLength_mSec);                
            end
        end
        
        % === 参数结构体返回函数 ===
        function simpleDelayParameters = getParameters(obj)
            simpleDelayParameters = obj.simpleDelayParameters;
        end
        
        % === 参数设置函数 ===
        % params = {delayTime_mSec, interpolate, delay_Samples}
        function setParameters(obj, params)
            obj.simpleDelayParameters = params;
            % 若延时时间比缓存区长，则重新分配缓存区空间
            if obj.simpleDelayParameters.delayTime_mSec >= obj.bufferLength_mSec
                obj.createDelayBuffer(obj.sampleRate, obj.simpleDelayParameters.delayTime_mSec); 
            end
            obj.simpleDelayParameters.delay_Samples = obj.simpleDelayParameters.delayTime_mSec * obj.samplesPerMSec;
            obj.delayBuffer.setInterpolate(obj.simpleDelayParameters.interpolate);
        end
        
        % === 模块处理数据函数 ===
        % 先读取数据，再写入数据
        function yn = processAudioSample(obj, xn)
            delay_Samples = obj.simpleDelayParameters.delay_Samples;
            % 延时为0，则直接返回输入数据
            if 0 == delay_Samples
                yn = xn;
                return;
            elseif delay_Samples < 1
                if abs(floor(delay_Samples)-delay_Samples) < 1e-8
                    offset = floor(delay_Samples);
                    if (0 == offset) 
                        yn = xn;
                        return;
                    else
                        yn = obj.delayBuffer.readBuffer(1);
                    end
                elseif "kNearest" == obj.simpleDelayParameters.interpolate
                    offset = floor(delayInSamples + 0.5);
                    if 0 == offset
                        yn = xn;
                        return;
                    else
                        yn = obj.delayBuffer.readBuffer(1);
                    end
                elseif "kLinear" == obj.simpleDelayParameters.interpolate
                    offset = floor(delayInSamples + 1);
                    fraction = offset - delayinSamples;
                    
                    sample1 = obj.delayBuffer.readBuffer(1);
                    sample2 = xn;

                    yn = circularbuffer_LinearInterpolation(sample1, sample2, fraction);
                else
                    offset = floorf(delayinSamples) + 1;
                    fraction = offset - delayinSamples;
                    sample0 = obj.delayBuffer.readBuffer(2);
                    sample1 = obj.delayBuffer.readBuffer(1);
                    sample2 = xn;
                    yn = circularbuffer_SquareInterpolation(sample0, sample1, sample2, fraction);
                end
            else
                yn = obj.delayBuffer.readBuffer(obj.simpleDelayParameters.delay_Samples);
            end
            
           
            obj.delayBuffer.writeBuffer(xn);
        end
        
        % === 延时缓存创建函数 ===
        function createDelayBuffer(obj, sampleRate_, bufferLength_mSec_)
            obj.bufferLength_mSec = bufferLength_mSec_;
            obj.sampleRate = sampleRate_;
            obj.samplesPerMSec = obj.sampleRate / 1000.0;
            obj.simpleDelayParameters.delay_Samples = obj.simpleDelayParameters.delayTime_mSec * obj.samplesPerMSec;
            obj.bufferLength = ceil(obj.bufferLength_mSec * obj.samplesPerMSec) + 1;
            obj.delayBuffer.createCircularBuffer(obj.bufferLength);
        end
        
        % === 读取延时数据函数 ===
        function output = readDelay(obj)
            output = obj.delayBuffer.readBuffer(obj.simpleDelayParameters.delay_Samples);
        end
        
        % === 读取一定时间前延时数据函数 ===
        function output = readDelayAtTime_mSec(obj, delay_mSec_)
            delay_Samples_ = delay_mSec_ * obj.samplesPerMSec;
            output = obj.delayBuffer.readBuffer(delay_Samples_);
        end
        
        % === 按照比例读取延时数据函数 ===
        function output = readDelayAtPercentage(obj, delayPercent)
            delay_Samples_ = (delayPercent / 100.0)*obj.simpleDelayParameters.delay_Samples;
            output = obj.delayBuffer.readBuffer(delay_Samples_);
        end
        
        % === 写延时数据函数 ===
        function writeDelay(obj, xn)
            obj.delayBuffer.writeBuffer(xn);
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