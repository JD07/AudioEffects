% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% 音频滤波器模块
% enum class filterAlgorithm {
%	kLPF1P, kLPF1, kHPF1, kLPF2, kHPF2, kBPF2, kBSF2, kButterLPF2, kButterHPF2, kButterBPF2,
%	kButterBSF2, kMMALPF2, kMMALPF2B, kLowShelf, kHiShelf, kNCQParaEQ, kCQParaEQ, kLWRLPF2, kLWRHPF2,
%	kAPF1, kAPF2, kResonA, kResonB, kMatchLP2A, kMatchLP2B, kMatchBP2A, kMatchBP2B,
%	kImpInvLP1, kImpInvLP2 } 
% 

classdef AudioFilter < handle
    % 利用常量类成员来模仿C++的枚举
    properties(Constant, Access = private)
        % enum filterCoeff { a0, a1, a2, b1, b2, c0, d0, numCoeffs }
        a0 = 1;
        a1 = 2;
        a2 = 3;
        b1 = 4;
        b2 = 5;
        c0 = 6;
        d0 = 7;
        numCoeffs = 7;
    end
    
    properties(Access = public)
        audioFilterParameters % 可调整参数结构体
        % --- algorithm % 滤波器类型
        % --- fc % 截止频率或中心频率
        % --- Q % 滤波器品质系数
        % --- boostCut_dB 滤波器增益，只在部分滤波器上使用
        biquad % 二阶滤波器类
        coeffArray % 系数矩阵
        sampleRate % 采样率
    end

    methods(Access = public)
        % === 构造函数 ===
        function obj = AudioFilter()
            obj.audioFilterParameters = struct();
            obj.audioFilterParameters.algorithm = "kLPF1";
            obj.audioFilterParameters.fc = 100.0;
            obj.audioFilterParameters.Q = 0.707;
            obj.audioFilterParameters.boostCut_dB = 0.0;

            obj.biquad = Biquad();
            obj.coeffArray = zeros(1, obj.numCoeffs);
            obj.sampleRate = 44100.0;
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            bqp = obj.biquad.getParameters();
            bqp.biquadCalcType = "kTransposeCanonical";
            obj.biquad.setParameters(bqp);
            obj.sampleRate = sampleRate_;
            obj.biquad.reset(sampleRate_);
        end
        
        % === 采样率设置函数 ===
        function setSampleRate(obj, sampleRate_)
            obj.sampleRate = sampleRate_;
            obj.calculateFilterCoeffs();
        end
        
        % === 获取参数 ===
        function audioFilterParameters = getParameters(obj)
            audioFilterParameters = obj.audioFilterParameters;
        end
        
        % === 参数设置函数 ===
        function setParameters(obj, parameters)
            if obj.audioFilterParameters.algorithm ~= parameters.algorithm || ...
               obj.audioFilterParameters.boostCut_dB ~= parameters.boostCut_dB || ...
               obj.audioFilterParameters.fc ~= parameters.fc || ...
               obj.audioFilterParameters.Q ~= parameters.Q
                obj.audioFilterParameters = parameters;
            else
                return;
            end
            
            if obj.audioFilterParameters.Q <= 0
                obj.audioFilterParameters.Q = 0.707;
            end
            obj.calculateFilterCoeffs();
        end
        
        % === 数据处理函数 ===
        function yn = processAudioSample(obj, xn)
            yn = obj.coeffArray(obj.d0) * xn + obj.coeffArray(obj.c0) * obj.biquad.processAudioSample(xn);
        end
        
    end
    

    
    methods(Access = protected)
        % === 根据滤波器类型，计算滤波器系数 ===
        function calculateFilterCoeffs(obj)
            % 原有系数清空
            obj.coeffArray(1 : obj.numCoeffs) = 0;
            % 赋予默认值
            obj.coeffArray(obj.a0) = 1.0;
            obj.coeffArray(obj.c0) = 1.0;
            obj.coeffArray(obj.d0) = 0.0;
            
            % 将可调参数结构体里的参数取出
            algorithm = obj.audioFilterParameters.algorithm;
            fc = obj.audioFilterParameters.fc;
            Q = obj.audioFilterParameters.Q;
            boostCut_dB = obj.audioFilterParameters.boostCut_dB;
            
            if "kImpInvLP1" == algorithm
                T = 1 / obj.sampleRate;
                omega = 2.0 * pi * fc;
                eT = exp(-T * omega);
                
                obj.coeffArray(obj.a0) = 1 - eT;
                obj.coeffArray(obj.a1) = 0.0;
                obj.coeffArray(obj.a2) = 0.0;
                obj.coeffArray(obj.b1) = -eT;
                obj.coeffArray(obj.b2) = 0.0;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kImpInvLP2" == algorithm
                alpha = 2.0 * pi * fc / obj.sampleRate;
                p_Re = -alpha / (2.0 * Q);
                zeta = 1.0 / (2.0 * Q);
                p_Im = alpha * sqrt(1 - zeta^2);
                c_Re = 0.0;    
                c_Im = alpha / (2.0 * sqrt(1 - zeta^2));
                eP_re = exp(p_Re);
                
                obj.coeffArray(obj.a0) = c_Re;
                obj.coeffArray(obj.a1) = -2.0 * (c_Re * cos(p_Im) + c_Im * sin(p_Im)) * exp(p_Re);
                obj.coeffArray(obj.a2) = 0.0;
                obj.coeffArray(obj.b1) = -2.0 * eP_re * cos(p_Im);
                obj.coeffArray(obj.b2) = eP_re^2;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kMatchLP2A" == algorithm
                theta_c = 2.0 * pi * fc / obj.sampleRate;
                
                q = 1.0 / (2.0 * Q);
                
                b_1 = 0.0;
                b_2 = exp(-2.0 * q * theta_c);
                
                if q <= 1.0
                    b_1 = -2.0 * exp(-q * theta_c) * cos(sqrt(1.0 - q^2) * theta_c);
                else
                    b_1 = -2.0 * exp(-q * theta_c) * cosh(sqrt(q^2 - 1.0) * theta_c);
                end
                
                B0 = (1.0 + b_1 + b_2) * (1.0 + b_1 + b_2);
                B1 = (1.0 - b_1 + b_2) * (1.0 - b_1 + b_2);
                B2 = -4.0 * b_2;
                
                phi_0 = 1.0 - sin(theta_c / 2.0) * sin(theta_c / 2.0);
                phi_1 = sin(theta_c / 2.0) * sin(theta_c / 2.0);
                phi_2 = 4.0* phi_0 * phi_1;
                
                R1 = (B0 * phi_0 + B1 * phi_1 + B2 * phi_2) * Q^2;
                A0 = B0;
                A1 = (R1 - A0 * phi_0) / phi_1;
                
                if A0 < 0.0
                    A0 = 0.0;
                end
                if A1 < 0.0
                    A1 = 0.0;
                end
                
                a_0 = 0.5 * (sqrt(A0) + sqrt(A1));
                a_1 = sqrt(A0) - a_0;
                a_2 = 0.0;
                
                obj.coeffArray(obj.a0) = a_0;
                obj.coeffArray(obj.a1) = a_1;
                obj.coeffArray(obj.a2) = a_2;
                obj.coeffArray(obj.b1) = b_1;
                obj.coeffArray(obj.b2) = b_2;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kMatchLP2B" == algorithm
                theta_c = 2.0 * pi * fc / obj.sampleRate;
                q = 1.0 / (2.0 * Q);
                
                b_1 = 0.0;
                b_2 = exp(-2.0 * q * theta_c);
                
                if q <= 1.0
                    b_1 = -2.0 * exp(-q * theta_c) * cos(sqrt(1.0 - q^2) * theta_c);
                else
                    b_1 = -2.0 * exp(-q * theta_c) * cosh(sqrt(q^2 - 1.0) * theta_c);
                end
                
                f0 = theta_c / pi;
                
                r0 = 1.0 + b_1 + b_2;
                denom = (1.0 - f0^2)^2 + f0^2 / Q^2;
                denom = sqrt(denom);
                r1 = (1.0 - b_1 + b_2)*f0^2 / denom;
                
                a_0 = (r0 + r1) / 2.0;
                a_1 = r0 - a_0;
                a_2 = 0.0;
                
                obj.coeffArray(obj.a0) = a_0;
                obj.coeffArray(obj.a1) = a_1;
                obj.coeffArray(obj.a2) = a_2;
                obj.coeffArray(obj.b1) = b_1;
                obj.coeffArray(obj.b2) = b_2;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kMatchBP2A" == algorithm
                theta_c = 2.0 * pi * fc / obj.sampleRate;
                q = 1.0 / (2.0 * Q);
                
                b_1 = 0.0;
                b_2 = exp(-2.0*q*theta_c);
                
                if q <= 1.0
                    b_1 = -2.0 * exp(-q * theta_c) * cos(sqrt(1 - q^2) * theta_c);
                else
                    b_1 = -2.0 * exp(-q * theta_c) * cosh(sqrt(q^2 - 1) * theta_c);
                end
                
                B0 = (1.0 + b_1 + b_2)^2;
                B1 = (1.0 - b_1 + b_2)^2;
                B2 = -4.0 * b_2;
                
                phi_0 = 1.0 - sin(theta_c / 2.0) * sin(theta_c / 2.0);
                phi_1 = sin(theta_c / 2.0) * sin(theta_c / 2.0);
                phi_2 = 4.0 * phi_0 * phi_1;
                
                R1 = B0 * phi_0 + B1 * phi_1 + B2 * phi_2;
                R2 = -B0 + B1 + 4.0 * (phi_0 - phi_1) * B2;
                
                A2 = (R1 - R2 * phi_1) / (4.0 * phi_1^2);
                A1 = R2 + 4.0 * (phi_1 - phi_0) * A2;
                
                a_1 = -0.5 * sqrt(A1);
                a_0 = 0.5 * (sqrt(A2 + a_1^2) - a_1);
                a_2 = -a_0 - a_1;
                
                obj.coeffArray(obj.a0) = a_0;
                obj.coeffArray(obj.a1) = a_1;
                obj.coeffArray(obj.a2) = a_2;
                obj.coeffArray(obj.b1) = b_1;
                obj.coeffArray(obj.b2) = b_2;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kMatchBP2B" == algorithm
                theta_c = 2.0 * pi * fc / obj.sampleRate;
                q = 1.0 / (2.0 * Q);
                
                b_1 = 0.0;
                b_2 = exp(-2.0 * q * theta_c);
                if q <= 1.0
                    b_1 = -2.0 * exp(-q * theta_c) * cos(sqrt(1.0 - q^2) * theta_c);
                else
                    b_1 = -2.0 * exp(-q * theta_c) * cosh(sqrt(q^2 - 1.0) * theta_c);
                end
                
                f0 = theta_c / pi;
                
                r0 = (1.0 + b_1 + b_2) / (pi * f0 * Q);
                denom = (1.0 - f0^2)^2 + f0^2 / Q^2;
                denom = sqrt(denom);
                
                r1 = (1.0 - b_1 + b_2) * (f0 / Q) / denom;
                
                a_1 = -r1 / 2.0;
                a_0 = (r0 - a_1) / 2.0;
                a_2 = -a_0 - a_1;
                
                obj.coeffArray(obj.a0) = a_0;
                obj.coeffArray(obj.a1) = a_1;
                obj.coeffArray(obj.a2) = a_2;
                obj.coeffArray(obj.b1) = b_1;
                obj.coeffArray(obj.b2) = b_2;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kLPF1P" == algorithm
                theta_c = 2.0 * pi * fc / obj.sampleRate;
                gamma = 2.0 - cos(theta_c);
                
                filter_b1 = sqrt(gamma^2 - 1) - gamma;
                filter_a0 = 1 + filter_b1;
                
                obj.coeffArray(obj.a0) = filter_a0;
                obj.coeffArray(obj.a1) = 0.0;
                obj.coeffArray(obj.a2) = 0.0;
                obj.coeffArray(obj.b1) = filter_b1;
                obj.coeffArray(obj.b2) = 0.0;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kLPF1" == algorithm
                theta_c = 2.0 * pi * fc / obj.sampleRate;
                gamma = cos(theta_c) / (1.0 + sin(theta_c));
                
                obj.coeffArray(obj.a0) = (1.0 - gamma) / 2.0;
                obj.coeffArray(obj.a1) = (1.0 - gamma) / 2.0;
                obj.coeffArray(obj.a2) = 0.0;
                obj.coeffArray(obj.b1) = -gamma;
                obj.coeffArray(obj.b2) = 0.0;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kHPF1" == algorithm
                theta_c = 2.0 * pi * fc / obj.sampleRate;
                gamma = cos(theta_c) / (1.0 + sin(theta_c));
                
                obj.coeffArray(obj.a0) = (1.0 + gamma) / 2.0;
                obj.coeffArray(obj.a1) = -(1.0 + gamma) / 2.0;
                obj.coeffArray(obj.a2) = 0.0;
                obj.coeffArray(obj.b1) = -gamma;
                obj.coeffArray(obj.b2) = 0.0;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kLPF2" == algorithm
                theta_c = 2.0 * pi * fc / obj.sampleRate;
                d = 1.0 / Q;
                betaNumerator = 1.0 - (d / 2.0 * sin(theta_c));
                betaDenominator = 1.0 + (d / 2.0 * sin(theta_c));
                
                beta = 0.5 * (betaNumerator / betaDenominator);
                gamma = (0.5 + beta) * cos(theta_c);
                alpha = (0.5 + beta - gamma) / 2.0;

                obj.coeffArray(obj.a0) = alpha;
                obj.coeffArray(obj.a1) = 2.0 * alpha;
                obj.coeffArray(obj.a2) = alpha;
                obj.coeffArray(obj.b1) = -2.0 * gamma;
                obj.coeffArray(obj.b2) = 2.0 * beta;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kHPF2" == algorithm
                theta_c = 2.0 * pi * fc / obj.sampleRate;
                d = 1.0 / Q;
                
                betaNumerator = 1.0 - (d / 2.0 * sin(theta_c));
                betaDenominator = 1.0 + (d / 2.0 * sin(theta_c));
                
                beta = 0.5 * betaNumerator / betaDenominator;
                gamma = (0.5 + beta) * cos(theta_c);
                alpha = (0.5 + beta + gamma) / 2.0;
                
                obj.coeffArray(obj.a0) = alpha;
                obj.coeffArray(obj.a1) = -2.0*alpha;
                obj.coeffArray(obj.a2) = alpha;
                obj.coeffArray(obj.b1) = -2.0*gamma;
                obj.coeffArray(obj.b2) = 2.0*beta;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kBPF2" == algorithm
                K = tan(pi * fc / obj.sampleRate);
                delta = K^2*Q + K + Q;
                
                obj.coeffArray(obj.a0) = K / delta;
                obj.coeffArray(obj.a1) = 0.0;
                obj.coeffArray(obj.a2) = -K / delta;
                obj.coeffArray(obj.b1) = 2.0*Q*(K^2 - 1) / delta;
                obj.coeffArray(obj.b2) = (K^2*Q - K + Q) / delta;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kBSF2" == algorithm
                K = tan(pi * fc / obj.sampleRate);
                delta = K^2 * Q + K + Q;
                
                obj.coeffArray(obj.a0) = Q * (1 + K^2) / delta;
                obj.coeffArray(obj.a1) = 2.0 * Q * (K^2 - 1) / delta;
                obj.coeffArray(obj.a2) = Q * (1 + K^2) / delta;
                obj.coeffArray(obj.b1) = 2.0*Q*(K^2 - 1) / delta;
                obj.coeffArray(obj.b2) = (K^2*Q - K + Q) / delta;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kButterLPF2" == algorithm
                theta_c = pi * fc / obj.sampleRate;
                C = 1.0 / tan(theta_c);
                
                obj.coeffArray(obj.a0) = 1.0 / (1.0 + sqrt(2)*C + C^2);
                obj.coeffArray(obj.a1) = 2.0 * obj.coeffArray(obj.a0);
                obj.coeffArray(obj.a2) = obj.coeffArray(obj.a0);
                obj.coeffArray(obj.b1) = 2.0 * obj.coeffArray(obj.a0) * (1.0 - C^2);
                obj.coeffArray(obj.b2) = obj.coeffArray(obj.a0) * (1.0 - sqrt(2)*C + C^2);
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kButterHPF2" == algorithm
                theta_c = pi * fc / obj.sampleRate;
                C = tan(theta_c);
                
                obj.coeffArray(obj.a0) = 1.0 / (1.0 + sqrt(2)*C + C^2);
                obj.coeffArray(obj.a1) = -2.0 * obj.coeffArray(obj.a0);
                obj.coeffArray(obj.a2) = obj.coeffArray(obj.a0);
                obj.coeffArray(obj.b1) = 2.0 * obj.coeffArray(obj.a0) * (C^2 - 1.0);
                obj.coeffArray(obj.b2) = obj.coeffArray(obj.a0) * (1.0 - sqrt(2)*C + C^2);
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kButterBPF2" == algorithm
                theta_c = 2.0 * pi * fc / obj.sampleRate;
                BW = fc / Q;
                delta_c = pi * BW / obj.sampleRate;
                delta_c = min(delta_c, 0.95 * pi / 2.0);
                
                C = 1.0 / tan(delta_c);
                D = 2.0 * cos(theta_c);
                
                obj.coeffArray(obj.a0) = 1.0 / (1.0 + C);
                obj.coeffArray(obj.a1) = 0.0;
                obj.coeffArray(obj.a2) = -obj.coeffArray(obj.a0);
                obj.coeffArray(obj.b1) = -obj.coeffArray(obj.a0) * C * D;
                obj.coeffArray(obj.b2) = obj.coeffArray(obj.a0) * (C - 1.0);
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kButterBSF2" == algorithm
                theta_c = 2.0 * pi * fc / obj.sampleRate;
                BW = fc / Q;
                delta_c = pi * BW / obj.sampleRate;
                delta_c = min(delta_c, 0.95 * pi / 2.0);

                C = tan(delta_c);
                D = 2.0 * cos(theta_c);
                
                obj.coeffArray(obj.a0) = 1.0 / (1.0 + C);
                obj.coeffArray(obj.a1) = -obj.coeffArray(obj.a0) * D;
                obj.coeffArray(obj.a2) = obj.coeffArray(obj.a0);
                obj.coeffArray(obj.b1) = -obj.coeffArray(obj.a0) * D;
                obj.coeffArray(obj.b2) = obj.coeffArray(obj.a0) * (1.0 - C);
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kMMALPF2" == algorithm || "kMMALPF2B" == algorithm
                theta_c = 2.0 * pi * fc / obj.sampleRate;
                resonance_dB = 0;
                
                if Q > 0.707
                    peak = Q^2 / sqrt(Q^2 - 0.25);
                    resonance_dB = 20.0 * log10(peak);
                end
                
                resonance = (cos(theta_c) + sin(theta_c) * sqrt(10^(resonance_dB/10) - 1))/...
                             (10^(resonance_dB/20.0) * sin(theta_c) + 1);
                g = 10 ^ (-resonance_dB / 40.0);
                         
                if "kMMALPF2B" == algorithm         
                    g = 1.0;
                end
                
                filter_b1 = -2.0 * resonance * cos(theta_c);
                filter_b2 = resonance^2;
                filter_a0 = g * (1 + filter_b1 + filter_b2);
                    
                obj.coeffArray(obj.a0) = filter_a0;
                obj.coeffArray(obj.a1) = 0.0;
                obj.coeffArray(obj.a2) = 0.0;
                obj.coeffArray(obj.b1) = filter_b1;
                obj.coeffArray(obj.b2) = filter_b2;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kLowShelf" == algorithm
                theta_c = 2.0 * pi * fc / obj.sampleRate;
                mu = 10 ^ (boostCut_dB / 20.0);
                
                beta = 4.0 / (1.0 + mu);
                delta = beta * tan(theta_c / 2.0);
                gamma = (1.0 - delta) / (1.0 + delta);
                
                obj.coeffArray(obj.a0) = (1.0 - gamma) / 2.0;
                obj.coeffArray(obj.a1) = (1.0 - gamma) / 2.0;
                obj.coeffArray(obj.a2) = 0.0;
                obj.coeffArray(obj.b1) = -gamma;
                obj.coeffArray(obj.b2) = 0.0;
                
                obj.coeffArray(obj.c0) = mu - 1.0;
                obj.coeffArray(obj.d0) = 1.0;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kHiShelf" == algorithm
                theta_c = 2.0 * pi * fc / obj.sampleRate;
                mu = 10 ^ (boostCut_dB / 20.0);
                
                beta = (1.0 + mu) / 4.0;
                delta = beta * tan(theta_c / 2.0);
                gamma = (1.0 - delta) / (1.0 + delta);
                
                obj.coeffArray(obj.a0) = (1.0 + gamma) / 2.0;
                obj.coeffArray(obj.a1) = -obj.coeffArray(obj.a0);
                obj.coeffArray(obj.a2) = 0.0;
                obj.coeffArray(obj.b1) = -gamma;
                obj.coeffArray(obj.b2) = 0.0;
                
                obj.coeffArray(obj.c0) = mu - 1.0;
                obj.coeffArray(obj.d0) = 1.0;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kCQParaEQ" == algorithm
                K = tan(pi * fc / obj.sampleRate);
                Vo = 10 ^ (boostCut_dB / 20.0);
                if boostCut_dB >=0
                    bBoost = true;
                else
                    bBoost = false;
                end
                
                d = 1.0 + (1.0 / Q) * K + K^2;
                e0 = 1.0 + (1.0 / (Vo * Q))*K + K^2;
                alpha = 1.0 + (Vo / Q)*K + K^2;
                beta = 2.0 * (K^2 - 1.0);
                gamma = 1.0 - (Vo / Q) * K + K^2;
                delta = 1.0 - (1.0 / Q)*K + K^2;
                eta = 1.0 - (1.0 / (Vo*Q))*K + K^2;
                
                if true == bBoost
                    obj.coeffArray(obj.a0) = alpha / d;
                    obj.coeffArray(obj.a1) = beta / d;
                    obj.coeffArray(obj.a2) = gamma / d;
                    obj.coeffArray(obj.b1) = beta / d;
                    obj.coeffArray(obj.b2) = delta / d;
                else
                    obj.coeffArray(obj.a0) = d / e0;
                    obj.coeffArray(obj.a1) = beta / e0;
                    obj.coeffArray(obj.a2) = delta / e0;
                    obj.coeffArray(obj.b1) = beta / e0;
                    obj.coeffArray(obj.b2) = eta / e0;
                end
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kNCQParaEQ" == algorithm
                theta_c = 2.0 * pi * fc / obj.sampleRate;
                mu = 10 ^ (boostCut_dB / 20.0);
                
                tanArg = theta_c / (2.0 * Q);
                tanArg = min(tanArg, 0.95 * pi / 2.0);
                
                zeta = 4.0 / (1.0 + mu);
                betaNumerator = 1.0 - zeta * tan(tanArg);
                betaDenominator = 1.0 + zeta * tan(tanArg);
                
                beta = 0.5 * betaNumerator / betaDenominator;
                gamma = (0.5 + beta) * cos(theta_c);
                alpha = 0.5 - beta;
                
                obj.coeffArray(obj.a0) = alpha;
                obj.coeffArray(obj.a1) = 0.0;
                obj.coeffArray(obj.a2) = -alpha;
                obj.coeffArray(obj.b1) = -2.0 * gamma;
                obj.coeffArray(obj.b2) = 2.0 * beta;
                
                obj.coeffArray(obj.c0) = mu - 1.0;
                obj.coeffArray(obj.d0) = 1.0;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kLWRLPF2" == algorithm
                omega_c = pi * fc;
                theta_c = pi * fc / obj.sampleRate;
                
                k = omega_c / tan(theta_c);
                denominator = k^2 + omega_c^2 + 2.0 * k * omega_c;
                b1_Num = -2.0 * k^2 + 2.0 * omega_c^2;
                b2_Num = -2.0 * k * omega_c + k^2 + omega_c^2;
                
                obj.coeffArray(obj.a0) = omega_c^2 / denominator;
                obj.coeffArray(obj.a1) = 2.0 * omega_c^2 / denominator;
                obj.coeffArray(obj.a2) = obj.coeffArray(obj.a0);
                obj.coeffArray(obj.b1) = b1_Num / denominator;
                obj.coeffArray(obj.b2) = b2_Num / denominator;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kLWRHPF2" == algorithm
                omega_c = pi * fc;
                theta_c = pi * fc / obj.sampleRate;
                
                k = omega_c / tan(theta_c);
                denominator = k^2 + omega_c^2 + 2.0*k*omega_c;
                b1_Num = -2.0*k^2 + 2.0*omega_c^2;
                b2_Num = -2.0*k*omega_c + k^2 + omega_c^2;
                
                obj.coeffArray(obj.a0) = k^2 / denominator;
                obj.coeffArray(obj.a1) = -2.0 * k^2 / denominator;
                obj.coeffArray(obj.a2) = obj.coeffArray(obj.a0);
                obj.coeffArray(obj.b1) = b1_Num / denominator;
                obj.coeffArray(obj.b2) = b2_Num / denominator;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kAPF1" == algorithm
                alphaNumerator = tan((pi*fc) / obj.sampleRate) - 1.0;
                alphaDenominator = tan((pi*fc) / obj.sampleRate) + 1.0;
                alpha = alphaNumerator / alphaDenominator;
                
                obj.coeffArray(obj.a0) = alpha;
                obj.coeffArray(obj.a1) = 1.0;
                obj.coeffArray(obj.a2) = 0.0;
                obj.coeffArray(obj.b1) = alpha;
                obj.coeffArray(obj.b2) = 0.0;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kAPF2" == algorithm
                theta_c = 2.0 * pi * fc / obj.sampleRate;
                BW = fc / Q;
                argTan = pi * BW / obj.sampleRate;
                if argTan >= 0.95 * pi / 2.0
                    argTan = 0.95 * pi / 2.0;
                end
                
                alphaNumerator = tan(argTan) - 1.0;
                alphaDenominator = tan(argTan) + 1.0;
                alpha = alphaNumerator / alphaDenominator;
                beta = -cos(theta_c);
                
                obj.coeffArray(obj.a0) = -alpha;
                obj.coeffArray(obj.a1) = beta * (1.0 - alpha);
                obj.coeffArray(obj.a2) = 1.0;
                obj.coeffArray(obj.b1) = beta * (1.0 - alpha);
                obj.coeffArray(obj.b2) = -alpha;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kResonA" == algorithm
                theta_c = 2.0 * pi *fc / obj.sampleRate;
                BW = fc / Q;
                filter_b2 = exp(-2.0 * pi * (BW / obj.sampleRate));
                filter_b1 = ((-4.0 * filter_b2) / (1.0 + filter_b2)) * cos(theta_c);
                filter_a0 = (1.0 - filter_b2) * sqrt(1.0 - filter_b1^2 / (4.0 * filter_b2));
                
                obj.coeffArray(obj.a0) = filter_a0;
                obj.coeffArray(obj.a1) = 0.0;
                obj.coeffArray(obj.a2) = 0.0;
                obj.coeffArray(obj.b1) = filter_b1;
                obj.coeffArray(obj.b2) = filter_b2;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            elseif "kResonB" == algorithm
                theta_c = 2.0 * pi * fc / obj.sampleRate;
                BW = fc / Q;
                filter_b2 = exp(-2.0 * pi * (BW / obj.sampleRate));
                filter_b1 = ((-4.0 * filter_b2) / (1.0 + filter_b2)) * cos(theta_c);
                filter_a0 = 1.0 - sqrt(filter_b2);
                
                obj.coeffArray(obj.a0) = filter_a0;
                obj.coeffArray(obj.a1) = 0.0;
                obj.coeffArray(obj.a2) = -filter_a0;
                obj.coeffArray(obj.b1) = filter_b1;
                obj.coeffArray(obj.b2) = filter_b2;
                obj.biquad.setCoefficients(obj.coeffArray);
                
            end
            
        end
    end
    
end