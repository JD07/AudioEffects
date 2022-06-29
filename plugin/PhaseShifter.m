% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% PhaseShifter

classdef PhaseShifter < handle
    properties(Constant, Access = private)
        PHASER_STAGES = 6;
        apf0_minF = 32.0;
        apf0_maxF = 1500.0;

        apf1_minF = 68.0;
        apf1_maxF = 3400.0;

        apf2_minF = 96.0;
        apf2_maxF = 4800.0;

        apf3_minF = 212.0;
        apf3_maxF = 10000.0;

        apf4_minF = 320.0;
        apf4_maxF = 16000.0;

        apf5_minF = 636.0;
        apf5_maxF = 20480.0;
    end
    
    properties(Access = public)
        parameters % 可调整参数结构体
        % --- lfoRate_Hz
        % --- lfoDepth_Pct
        % --- intensity_Pct
        % --- quadPhaseLFO
        apf % AudioFilter矩阵
        lfo % LFO
    end
    
    methods(Access = public)
        % === 构造函数 ===
        function obj = PhaseShifter()
            for i = 1:obj.PHASER_STAGES
                obj.apf(i) = AudioFilter();
            end
            obj.lfo = LFO();
            
            lfoparams = obj.lfo.getParameters();
            lfoparams.waveform = "kTriangle";
            lfo.setParameters(lfoparams);
            
            params = obj.apf(0).getParameters();
            params.algorithm = "kAPF1";
            
            for i = 1:obj.PHASER_STAGES
                obj.apf(i).setParameters(params);
            end
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            obj.lfo.reset(sampleRate_);
            
            for i = 1:obj.PHASER_STAGES
                obj.apf(i).reset(sampleRate_);
            end
            
        end
        
        % === 数据处理函数 ===
        function output = processAudioSample(obj, xn)
            lfoData = obj.lfo.renderAudioOutput();
            
            lfoValue = lfoData.normalOutput;
            if (obj.parameters.quadPhaseLFO)
                lfoValue = lfoData.quadPhaseOutput_pos;
            end
            
            depth = obj.parameters.lfoDepth_Pct / 100.0;
            modulatorValue = lfoValue * depth;
            
            params = obj.apf(1).getParameters();
            params.fc = doBipolarModulation(modulatorValue, obj.apf0_minF, obj.apf0_maxF);
            obj.apf(1).setParameters(params);
            
            params = obj.apf(2).getParameters();
            params.fc = doBipolarModulation(modulatorValue, obj.apf1_minF, obj.apf1_maxF);
            obj.apf(2).setParameters(params);
            
            params = obj.apf(3).getParameters();
            params.fc = doBipolarModulation(modulatorValue, obj.apf2_minF, obj.apf2_maxF);
            obj.apf(3).setParameters(params);
            
            params = obj.apf(4).getParameters();
            params.fc = doBipolarModulation(modulatorValue, obj.apf3_minF, obj.apf3_maxF);
            obj.apf(4).setParameters(params);
            
            params = obj.apf(5).getParameters();
            params.fc = doBipolarModulation(modulatorValue, obj.apf4_minF, obj.apf4_maxF);
            obj.apf(5).setParameters(params);
            
            params = obj.apf(6).getParameters();
            params.fc = doBipolarModulation(modulatorValue, obj.apf5_minF, obj.apf5_maxF);
            obj.apf(6).setParameters(params);
            
            gamma1 = obj.apf(6).getG_value();
            gamma2 = obj.apf(5).getG_value()*gamma1;
            gamma3 = obj.apf(4).getG_value()*gamma2;
            gamma4 = obj.apf(3).getG_value()*gamma3;
            gamma5 = obj.apf(2).getG_value()*gamma4;
            gamma6 = obj.apf(1).getG_value()*gamma5;
            
            K = obj.parameters.intensity_Pct / 100.0;
            alpha0 = 1.0 / (1.0 + K * gamma6);
            
            Sn = gamma5*obj.apf(1).getS_value() + gamma4*obj.apf(2).getS_value() + gamma3*obj.apf(3).getS_value() +...
                gamma2*obj.apf(4).getS_value() + gamma1*obj.apf(5).getS_value() + obj.apf(6).getS_value();
            
            u = alpha0 * (xn + K*Sn);
            
            APF1 = obj.apf(1).processAudioSample(u);
            APF2 = obj.apf(2).processAudioSample(APF1);
            APF3 = obj.apf(3).processAudioSample(APF2);
            APF4 = obj.apf(4).processAudioSample(APF3);
            APF5 = obj.apf(5).processAudioSample(APF4);
            APF6 = obj.apf(6).processAudioSample(APF5);
            
            output = 0.125*xn + 1.25*APF6;
        end
        
        % === 获取参数 ===
        function parameters = getParameters(obj)
            parameters = obj.parameters;
        end
        
        % === 设置参数 ===
        function setParameters(obj, params)
            if params.lfoRate_Hz ~= obj.parameters.lfoRate_Hz
                lfoparams = obj.lfo.getParameters();
                lfoparams.frequency_Hz = params.lfoRate_Hz;
                obj.lfo.setParameters(lfoparams);
            end
            
            obj.parameters = params;
            
        end

    end
    
end