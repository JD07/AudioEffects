% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% 梳状滤波器模块
% 使用方法
% （1）调用构造函数初始化类的实例
% （2）调用reset方法设置采样率
% （3）调用getParameters方法，返回可调整参数结构体
% （4）设置可调整参数结构体，并利用setParameters方法设置参数

classdef CombFilter < handle
    properties(Access = public)
        combFilterParameters % 可调参数结构体
        % --- delayTime_mSec % 延迟时间（单位：ms）
        % --- RT60Time_mSec % 混响时间T60
        % --- enableLPF % 使能LPF标志
        % --- lpf_g % LPF增益值
        % --- interpolate % 插值标志
        sampleRate 
        comb_g % 梳状滤波器系数
        bufferLength_mSec % 缓存区长度
        
        lpf_g % LPF增益值
        lpf_state % LPF状态值
        
        delay % SimpleDelay类
    end
    
    methods(Access = public)
        % === 构造函数 ===
        function obj = CombFilter()
            obj.combFilterParameters = struct();
            obj.combFilterParameters.delayTime_mSec = 0.0;
            obj.combFilterParameters.RT60Time_mSec = 0.0;
            obj.combFilterParameters.enableLPF = false;
            obj.combFilterParameters.lpf_g = 0.0;
            obj.combFilterParameters.interpolate = "kNearest";

            obj.sampleRate = 0.0;
            obj.comb_g = 0.0;
            obj.bufferLength_mSec = 0.0;
            
            obj.lpf_g = 0.0;
            obj.lpf_state = 0.0;
            
            obj.delay = SimpleDelay();
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            obj.lpf_state = 0.0;
            obj.delay.reset(sampleRate_);
            obj.createDelayBuffer(sampleRate_, obj.bufferLength_mSec);
        end
        
        % === 数据处理函数 ===
        function yn = processAudioSample(obj, xn)
            % 延时为0，直接返回
            delayParams = obj.delay.getParameters();
            if (delayParams.delay_Samples == 0)
                yn = xn;
                return
            end
            
            yn = obj.delay.readDelay();
            
            if obj.combFilterParameters.enableLPF
                g2 = obj.lpf_g * (1 - obj.comb_g);
                filteredSignal = yn + g2 * obj.lpf_state;
                input = xn + obj.comb_g * filteredSignal;
                obj.lpf_state = filteredSignal;
            else
                input = xn + obj.comb_g * yn;
            end
            
            obj.delay.writeDelay(input);
        end
        
        % === 获取参数 ===
        function combFilterParameters = getParameters(obj)
            combFilterParameters = obj.combFilterParameters;
        end
        
        % === 参数设置函数 ===
        function setParameters(obj, params)
            obj.combFilterParameters = params;
            
            % 更新延迟模块参数
            delayParams = obj.delay.getParameters();
            delayParams.delayTime_mSec = obj.combFilterParameters.delayTime_mSec;
            delayParams.interpolate = obj.combFilterParameters.interpolate;
            obj.delay.setParameters(delayParams);
            
            % 计算T60参数
            delayParams = obj.delay.getParameters();
            exponent = -3.0 * delayParams.delay_Samples * (1.0 / obj.sampleRate);
            rt60_mSec = obj.combFilterParameters.RT60Time_mSec / 1000.0;
            obj.comb_g = 10 ^ (exponent / rt60_mSec);
            
            obj.lpf_g = obj.combFilterParameters.lpf_g;
        end
        
        function createDelayBuffer(obj, sampleRate_, delay_mSec)
            obj.sampleRate = sampleRate_;
            obj.bufferLength_mSec = delay_mSec;
            % 创建缓存区
            obj.delay.createDelayBuffer(sampleRate_, delay_mSec);
        end
        
    end
end