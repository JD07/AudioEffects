% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% 延迟全通滤波器模块

classdef DelayAPF < handle
    properties(Access = public)
        delayAPFParameters % 可调制参数结构体
        % --- delayTime_mSec % APF延迟时间
        % --- apf_g % APF增益系数
        % --- enableLPF % LPF使能标志
        % --- lpf_g % LPF增益系数
        % --- interpolate % 插值标志
        % --- enableLFO % LFO使能标志
        % --- lfoRate_Hz % LFO频率
        % --- lfoDepth % LFO范围
        % --- lfoMaxModulation_mSec % LFO最大调制时间
        sampleRate 
        bufferLength_mSec
        delay % SimpleDelay类
        modLFO % LFO类
        lpf_state % lpf缓存
    end
    
    methods(Access = public)
        % === 构造函数 ===
        function obj = DelayAPF()
            obj.delayAPFParameters = struct();
            obj.delayAPFParameters.delayTime_mSec = 0.0;
            obj.delayAPFParameters.apf_g = 0.0;
            obj.delayAPFParameters.enableLPF = false;
            obj.delayAPFParameters.lpf_g = 0.0;
            obj.delayAPFParameters.interpolate = "kNearest"; % 插值类型 kNearest, kLinear, kSquare, kCubic;
            obj.delayAPFParameters.enableLFO = false;
            obj.delayAPFParameters.lfoRate_Hz = 0.0;
            obj.delayAPFParameters.lfoDepth = 0.0;
            obj.delayAPFParameters.lfoMaxModulation_mSec = 0.0;

            obj.sampleRate = 0.0;
            obj.bufferLength_mSec = 0.0;
            obj.delay = SimpleDelay();
            obj.modLFO = LFO();
            obj.lpf_state = 0.0;
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            obj.modLFO.reset(sampleRate_);
            obj.lpf_state = 0.0;
            obj.createDelayBuffer(sampleRate_, obj.bufferLength_mSec);
        end
        
        % === 数据处理函数 ===
        function yn = processAudioSample(obj, xn)
            delayParams = obj.delay.getParameters();
            if 0 == delayParams.delay_Samples
                yn = xn;
                return;
            end
            
            wnD = 0.0;
            apf_g = obj.delayAPFParameters.apf_g;
            lpf_g = obj.delayAPFParameters.lpf_g;
            lfoDepth = obj.delayAPFParameters.lfoDepth;
            
            if obj.delayAPFParameters.enableLFO
                lfoOutput = obj.modLFO.renderAudioOutput();
                maxDelay = delayParams.delayTime_mSec;
                minDelay = maxDelay - obj.delayAPFParameters.lfoMaxModulation_mSec;
                minDelay = max(0.0, minDelay);
            
                modDelay_mSec = doUnipolarModulationFromMax(bipolarToUnipolar(lfoDepth * lfoOutput.normalOutput),...
                                                            minDelay, maxDelay);
                wnD = obj.delay.readDelayAtTime_mSec(modDelay_mSec);
            else
                wnD = obj.delay.readDelay();
            end
            
            if obj.delayAPFParameters.enableLPF
                wnD = wnD*(1.0 - lpf_g) + lpf_g * obj.lpf_state;
                obj.lpf_state = wnD;
            end
            
            wn = xn + apf_g * wnD;
            
            yn = -apf_g * wn + wnD;
            
            obj.delay.writeDelay(wn);
            
        end
        
        % === 获取参数 ===
        function delayAPFParameters = getParameters(obj)
            delayAPFParameters = obj.delayAPFParameters;
        end
        
        % === 参数设置函数 ===
        function setParameters(obj, params)
            obj.delayAPFParameters = params;
            
            % 更新延时
            delayParams = obj.delay.getParameters();
            delayParams.delayTime_mSec = obj.delayAPFParameters.delayTime_mSec;
            obj.delay.setParameters(delayParams);
            
            % 更新LFO
            lfoParams = obj.modLFO.getParameters();
            lfoParams.frequency_Hz = obj.delayAPFParameters.lfoRate_Hz;
            lfoParams.waveform = 'kSin';
            obj.modLFO.setParameters(lfoParams);
            
        end
        
        % === 缓存区创建函数 ===
        
        function createDelayBuffer(obj, sampleRate_, delay_mSec)
            obj.sampleRate = sampleRate_;
            obj.bufferLength_mSec = delay_mSec;
            
            obj.delay.createDelayBuffer(sampleRate_, delay_mSec);
        end
        
    end
end

function output = doUnipolarModulationFromMax(unipolarModulatorValue, minValue, maxValue)
    unipolarModulatorValue = max(0.0, unipolarModulatorValue);
    unipolarModulatorValue = min(1.0, unipolarModulatorValue);
    
    output = unipolarModulatorValue * (maxValue - minValue) + minValue;
end

function output = bipolarToUnipolar(value)
    output = 0.5 * value + 0.5;
end