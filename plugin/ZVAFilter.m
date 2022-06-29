% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% ZVAFilter

classdef ZVAFilter < handle
    properties(Access = public)
        zvaFilterParameters % 可调参数结构体
        % --- filterAlgorithm % 算法类型 kLPF1 kHPF1 kAPF1 kSVF_LP kSVF_HP kSVF_BP kSVF_BS
        % --- fc
        % --- Q
        % --- filterOutputGain_dB
        % --- enableGainComp % 增益补偿使能标志位
        % --- matchAnalogNyquistLPF
        % --- selfOscillate
        % --- enableNLP
        sampleRate
        integrator_z
        alpha0
        alpha
        rho
        beta
        analogMatchSigma
    end
    
    methods(Access = public)
        % === 构造函数 ===
        function obj = ZVAFilter()
            field1 = "filterAlgorithm";
            value1 = "kSVF_LP";
            field2 = "fc";
            value2 = 1000.0;
            field3 = "Q";
            value3 = 0.707;
            field4 = "filterOutputGain_dB";
            value4 = 0.0;
            field5 = "enableGainComp";
            value5 = false;
            field6 = "matchAnalogNyquistLPF";
            value6 = false;
            field7 = "selfOscillate";
            value7 = false;
            field8 = "enableNLP";
            value8 = false;
            obj.zvaFilterParameters = struct(field1,value1,field2,value2,field3,value3,...
                 field4,value4,field5,value5,field6,value6,field7,value7,field8,value8);
            
            obj.sampleRate = 441000.0;
            obj.integrator_z = zeros(1,2);
            obj.alpha0 = 0.0;
            obj.alpha = 0.0;
            obj.rho = 0.0;
            
            obj.beta = 0.0;
            
            obj.analogMatchSigma = 0.0;
             
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            obj.sampleRate = sampleRate_;
            obj.integrator_z(1) = 0.0;
            obj.integrator_z(2) = 0.0;
        end
        
        % === 获取参数 ===
        function zvaFilterParameters = getParameters(obj)
            zvaFilterParameters = obj.zvaFilterParameters;
        end
        
        % === 参数设置 ===
        function setParameters(obj, params)
            if params.fc ~= obj.zvaFilterParameters.fc ||...
               params.Q ~= obj.zvaFilterParameters.Q ||...
               params.selfOscillate ~= obj.zvaFilterParameters.selfOscillate||...
               params.matchAnalogNyquistLPF ~= obj.zvaFilterParameters.matchAnalogNyquistLPF
                
                obj.zvaFilterParameters = params;
                obj.calculateFilterCoeffs();
            else
                obj.zvaFilterParameters = params;
            end
        end
        
        % === 数据处理函数 ===
        function yn = processAudioSample(obj, xn)
            filterAlgorithm = obj.zvaFilterParameters.filterAlgorithm;
            matchAnalogNyquistLPF = obj.zvaFilterParameters.matchAnalogNyquistLPF;
            
            if obj.zvaFilterParameters.enableGainComp
                peak_dB = obj.dBPeakGainFor_Q(obj.zvaFilterParameters.Q);
                if peak_dB > 0.0
                    halfPeak_dBGain = dB2Raw(-peak_dB / 2.0);
                    xn = xn * halfPeak_dBGain;
                end
            end
            
            if "kLPF1" == filterAlgorithm||...
               "kHPF1" == filterAlgorithm||...
               "kAPF1" == filterAlgorithm
                
                vn = (xn - obj.integrator_z(1)) * obj.alpha;
                lpf = (xn - obj.integrator_z(1)) * obj.alpha + obj.integrator_z(1);
                
                obj.integrator_z(1) = vn + lpf;
                hpf = xn - lpf;
                apf = lpf - hpf;
                
                if "kLPF1" == filterAlgorithm
                    if matchAnalogNyquistLPF
                        yn = lpf + alpha * hpf;
                        return;
                    else
                        yn = lpf;
                        return;
                    end
                elseif "kHPF1" == filterAlgorithm
                    yn = hpf;
                    return;
                elseif "kAPF1" == filterAlgorithm
                    yn = apf;
                    return;
                end
                
                yn = xn;
                return;
            end
            
            hpf = obj.alpha0 * (xn - obj.rho * obj.integrator_z(1) - obj.integrator_z(2));
            bpf = obj.alpha * hpf + obj.integrator_z(1);
            
            if obj.zvaFilterParameters.enableNLP
                bpf = softClipWaveShaper(bpf, 1.0);
            end
            
            lpf = obj.alpha * bpf + obj.integrator_z(2);
            bsf = hpf + lpf;
            sn = obj.integrator_z(1);
            
            obj.integrator_z(1) = obj.alpha * hpf + bpf;
            obj.integrator_z(2) = obj.alpha * bpf + lpf;
            
            filterOutputGain = 10^(obj.zvaFilterParameters.filterOutputGain_dB / 20.0);
            
            if "kSVF_LP" == filterAlgorithm
                if matchAnalogNyquistLPF
                    lpf = lpf + obj.analogMatchSigma * sn;
                end
                yn = filterOutputGain * lpf;
                return;
            elseif "kSVF_HP" == filterAlgorithm
                yn = filterOutputGain * hpf;
                return;
            elseif "kSVF_BP" == filterAlgorithm
                yn = filterOutputGain * bpf;
                return;
            elseif "kSVF_BS" == filterAlgorithm
                yn = filterOutputGain*bsf;
                return;
            end
            
            yn = filterOutputGain*lpf;
            return;
        end
        
        
        % === 计算滤波器系数 ===
        function calculateFilterCoeffs(obj)
            fc = obj.zvaFilterParameters.fc;
            Q = obj.zvaFilterParameters.Q;
            filterAlgorithm = obj.zvaFilterParameters.filterAlgorithm;
            
            wd = 2 * pi * fc;
            T = 1.0 / obj.sampleRate;
            wa = 2.0 / T * tan(wd * T / 2.0);
            g = wa * T / 2.0;
            
            if "kLPF1" == filterAlgorithm || "kHPF1" == filterAlgorithm || "kAPF1" == filterAlgorithm
                obj.alpha = g / (1.0 + g);
            else
                if obj.zvaFilterParameters
                    R = 0.0;
                else
                    R = 1.0 / (2.0 * Q);
                end
                
                obj.alpha0 = 1.0 / (1.0 + 2.0*R*g + g^2);
                obj.alpha = g;
                obj.rho = 2.0 * R + g;
                
                
            end
        end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    end
end