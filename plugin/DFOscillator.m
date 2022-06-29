% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% DFOscillator

classdef DFOscillator < handle
    properties(Constant, Access = private)
        % enum DFOscillatorCoeffs { df_b1, df_b2, numDFOCoeffs }
        df_b1 = 1;
        df_b2 = 2;
        numDFOCoeffs = 2;
        % enum DFOscillatorStates { df_yz1, df_yz2, numDFOStates }
        df_yz1 = 1;
        df_yz2 = 2;
        numDFOStates = 2;
    end
    
    properties(Access = public)
        parameters % 可调参数结构体
        % --- waveform % kTriangle, kSin, kSaw
        % --- frequency_Hz
        stateArray
        coeffArray
        sampleRate
    end
    
    methods
        % === 构造函数 ===
        function obj = DFOscillator()
            obj.parameters = struct();
            obj.parameters.waveform = "kTriangle";
            obj.parameters.frequency_Hz = 0.0;
            
            obj.stateArray = [0.0, 0.0];
            obj.coeffArray = [0.0, 0.0];
            
            obj.sampleRate = 0.0;
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            obj.sampleRate = sampleRate_;
            obj.stateArray(:) = 0;
            obj.updateDFO();
        end
        
        % === 获取参数 ===
        function parameters = getParameters(obj)
            parameters = obj.parameters;
        end
        
        % === 设置参数 ===
        function setParameters(obj, params)
            if obj.parameters.frequency_Hz ~= params.frequency_Hz
                obj.parameters = params;
                obj.updateDFO();
            end
        end
        
        % === 获取振荡器输出 ===
        function output = renderAudioOutput(obj)
            output = struct();
            
            % 差分方程计算
            output.normalOutput = -obj.coeffArray(obj.df_b1)*obj.stateArray(obj.df_yz1)...
                                  -obj.coeffArray(obj.df_b2)*obj.stateArray(obj.df_yz2);
            output.invertedOutput = -output.normalOutput;
            
            % 状态更新
            obj.stateArray(obj.df_yz2) = obj.stateArray(obj.df_yz1);
            obj.stateArray(obj.df_yz1) = output.normalOutput;
        end
        
        % === DFO更新 ===
        function updateDFO(obj)
            wT = (2*pi*obj.parameters.frequency_Hz) / obj.sampleRate;
            
            obj.coeffArray(obj.df_b1) = -2.0*cos(wT);
            obj.coeffArray(obj.df_b2) = 1.0;
            
            wnT1 = asin(obj.stateArray(obj.df_yz1));
            n = wnT1 / wT;
            
            if obj.stateArray(obj.df_yz1) > obj.stateArray(obj.df_yz2)
                n = n - 1;
            else
                n = n + 1;
            end
            
            obj.stateArray(obj.df_yz2) = sin(n * wT);
        end
        
    end
end