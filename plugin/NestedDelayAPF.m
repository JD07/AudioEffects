% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% 书本的NestedDelayAPF模块利用继承来实现嵌套
% 但是matlab的继承和C++有些区别，且继承不方便改写为C，所以在这里将代码进行修改

classdef NestedDelayAPF < handle
    properties(Access = public)
        nestedAPFParameters % 可调整参数结构体
        % --- outerAPFdelayTime_mSec % 外部apf延时时间
        % --- innerAPFdelayTime_mSec % 内部apf延时时间
        % --- outerAPF_g % 外部apf增益系数
        % --- innerAPF_g % 内部apf增益系数
        % --- enableLPF % LPF使能标志
        % --- lpf_g % LPF增益系数
        % --- interpolate % 插值标志
        % --- enableLFO % LFO使能标志
        % --- lfoRate_Hz % LFO震荡频率
        % --- lfoDepth % LFO震荡范围
        % --- lfoMaxModulation_mSec % 最大调制时间
        
        sampleRate 
        bufferLength_mSec
        delay % SimpleDelay类
        modLFO % LFO类
        lpf_state % lpf缓存
        
        innerDelayAPF % DelayAPF类
    end
   
    methods(Access = public)
        % === 构造函数 ===
        function obj = NestedDelayAPF()
            obj.nestedAPFParameters = struct();
            obj.nestedAPFParameters.outerAPFdelayTime_mSec = 0.0;
            obj.nestedAPFParameters.innerAPFdelayTime_mSec = 0.0;
            obj.nestedAPFParameters.outerAPF_g = 0.0;
            obj.nestedAPFParameters.innerAPF_g = 0.0;
            obj.nestedAPFParameters.enableLPF = false;
            obj.nestedAPFParameters.lpf_g = 0.0;
            obj.nestedAPFParameters.interpolate = "kNearest"; % 插值类型 kNearest, kLinear, kSquare, kCubic;
            obj.nestedAPFParameters.enableLFO = false;
            obj.nestedAPFParameters.lfoRate_Hz = 0.0;
            obj.nestedAPFParameters.lfoDepth = 0.0;
            obj.nestedAPFParameters.lfoMaxModulation_mSec = 0.0;
            
            obj.sampleRate = 0.0;
            obj.bufferLength_mSec = 0.0;
            obj.delay = SimpleDelay();
            obj.modLFO = LFO();
            obj.lpf_state = 0.0;
            
            obj.innerDelayAPF = DelayAPF();

        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            obj.sampleRate = sampleRate_;
            obj.modLFO.reset(sampleRate_);
            obj.lpf_state = 0.0;
            obj.innerDelayAPF.reset(sampleRate_);
        end
        
        % === 数据处理函数 ===
        function yn = processAudioSample(obj, xn)
            % --- delay line output ---
            wnD = 0.0;
            
            if 0 == obj.nestedAPFParameters.outerAPFdelayTime_mSec
                yn = obj.innerDelayAPF.processAudioSample(xn);
                return;
            end
            
            apf_g = obj.nestedAPFParameters.outerAPF_g;
            lpf_g = obj.nestedAPFParameters.lpf_g;
            lfoDepth = obj.nestedAPFParameters.lfoDepth;
            
            if obj.nestedAPFParameters.enableLFO
                lfoOutput = obj.modLFO.renderAudioOutput();
                maxDelay = obj.nestedAPFParameters.outerAPFdelayTime_mSec;
                minDelay = maxDelay - obj.nestedAPFParameters.lfoMaxModulation_mSec;
                minDelay = max(0.0, minDelay);
                
                modDelay_mSec = doUnipolarModulationFromMax(...
                    bipolarToUnipolar(lfoDepth*lfoOutput.normalOutput),minDelay, maxDelay);
                
                wnD = obj.delay.readDelayAtTime_mSec(modDelay_mSec);
            else
                wnD = obj.delay.readDelay();
            end
            
            if obj.nestedAPFParameters.enableLPF
                wnD = wnD * (1.0 - lpf_g) + lpf_g * obj.lpf_state;
                obj.lpf_state = wnD;
            end
            
            wn = xn + apf_g * wnD;
            
            ynInner = obj.innerDelayAPF.processAudioSample(wn);
            
            yn = -apf_g * wn + wnD;
            
            obj.delay.writeDelay(ynInner);
            
        end
        
        % === 获取参数 ===
        function nestedAPFParameters = getParameters(obj)
            nestedAPFParameters = obj.nestedAPFParameters;
        end
        
        % === 参数设置函数 ===
        function setParameters(obj, params)
            obj.nestedAPFParameters = params;
            
            % 更新延时
            delayParams = obj.delay.getParameters();
            delayParams.delayTime_mSec = obj.nestedAPFParameters.outerAPFdelayTime_mSec;
            obj.delay.setParameters(delayParams);
            
            % 更新LFO
            lfoParams = obj.modLFO.getParameters();
            lfoParams.frequency_Hz = obj.nestedAPFParameters.lfoRate_Hz;
            lfoParams.waveform = 'kSin';
            obj.modLFO.setParameters(lfoParams);
            
            % 更新DelayAPF
            delayAPFParams = obj.innerDelayAPF.getParameters();
            delayAPFParams.delayTime_mSec = obj.nestedAPFParameters.innerAPFdelayTime_mSec;
            delayAPFParams.apf_g = obj.nestedAPFParameters.innerAPF_g;
            delayAPFParams.enableLPF = false;
            delayAPFParams.lpf_g = 0.0;
            delayAPFParams.interpolate = "kNearest"; % 插值类型 kNearest, kLinear, kSquare, kCubic;
            delayAPFParams.enableLFO = false;
            delayAPFParams.lfoRate_Hz = 0.0;
            delayAPFParams.lfoDepth = 0.0;
            delayAPFParams.lfoMaxModulation_mSec = 0.0;
            obj.innerDelayAPF.setParameters(delayAPFParams);
            
        end
        
        % === 缓存区创建函数 ===
        function createDelayBuffers(obj, sampleRate_, delay_mSec_out, delay_mSec_in)
            obj.sampleRate = sampleRate_;
            obj.bufferLength_mSec = delay_mSec_out;
            
            obj.delay.createDelayBuffer(sampleRate_, delay_mSec_out);
            obj.innerDelayAPF.createDelayBuffer(sampleRate_, delay_mSec_in);
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