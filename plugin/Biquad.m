% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% 二阶滤波器模块
% enum filterCoeff { a0, a1, a2, b1, b2, c0, d0, numCoeffs };
% enum stateReg { x_z1, x_z2, y_z1, y_z2, numStates };


classdef Biquad < handle
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
        % enum stateReg { x_z1, x_z2, y_z1, y_z2, numStates }
        x_z1 = 1;
        x_z2 = 2;
        y_z1 = 3;
        y_z2 = 4;
        numStates = 4;
    end
    
    properties(Access = public)
        parameters
        % --- biquadCalcType % 滤波器结构 : kDirect kCanonical kTransposeDirect kTransposeCanonical
        coeffArray
        stateArray
        storageComponent
    end
    
    methods(Access = public)
        % === 构造函数 ===
        function obj = Biquad()
            obj.parameters = struct();
            obj.parameters.biquadCalcType = "kDirect";
            obj.coeffArray = zeros(1, obj.numCoeffs);
            obj.stateArray = zeros(1, obj.numStates);
            obj.storageComponent = 0.0;
        end
        
        % === 重置函数 ===
        function reset(obj, ~)
            obj.stateArray(1 : obj.numStates) = 0;
        end
        
        % === 参数获取 ===
        function parameters = getParameters(obj)
            parameters = obj.parameters;
        end
        
        % === 设置参数 ===
        function setParameters(obj, parameters_)
            obj.parameters = parameters_;
        end
        
        % === 设置系数 ===
        % 输入的系数是矩阵，源代码是使用指针和内存拷贝实现的设置，matlab可以直接赋值
        function setCoefficients(obj, coeffs)
            obj.coeffArray = coeffs;
        end
        
        % === 获取系数矩阵 ===
        function coeffArray = getCoefficients(obj)
            coeffArray = obj.coeffArray;
        end
        
        % === 获取状态矩阵 ===
        function stateArray = getStateArray(obj)
            stateArray = obj.stateArray;
        end
        
        % === 为Harma滤波器获取G(gain)值 ===
        function G = getG_value(obj)
            G = obj.coeffArray(obj.a0);
        end
        
        % === 为Harma滤波器获取S(storage)值 ===
        function storageComponent = getS_value(obj)
            if "kDirect" == obj.parameters.biquadCalcType
                storageComponent = obj.coeffArray(obj.a1) * obj.stateArray(obj.x_z1) +...
                                   obj.coeffArray(obj.a2) * obj.stateArray(obj.x_z2) -...
                                   obj.coeffArray(obj.b1) * obj.stateArray(obj.y_z1) -...
                                   obj.coeffArray(obj.b2) * obj.stateArray(obj.y_z2);
            elseif "kTransposeCanonical" == obj.parameters.biquadCalcType
                storageComponent = obj.stateArray(obj.x_z1);
            end
        end
        
        % === 处理函数 ===
        function yn = processAudioSample(obj, xn)
            % 直接型
            % y[n] = a0*x[n] + a1*x[n-1] + a2*x[n-2] - b1*y[n-1] - b2*y[n-2]
            if "kDirect" == obj.parameters.biquadCalcType
                yn = obj.coeffArray(obj.a0) * xn + ...
                     obj.coeffArray(obj.a1) * obj.stateArray(obj.x_z1) + ...
                     obj.coeffArray(obj.a2) * obj.stateArray(obj.x_z2) - ...
                     obj.coeffArray(obj.b1) * obj.stateArray(obj.y_z1) - ...
                     obj.coeffArray(obj.b2) * obj.stateArray(obj.y_z2);
                 % 状态更新
                 obj.stateArray(obj.x_z2) = obj.stateArray(obj.x_z1);
                 obj.stateArray(obj.x_z1) = xn;
                 
                 obj.stateArray(obj.y_z2) = obj.stateArray(obj.y_z1);
                 obj.stateArray(obj.y_z1) = yn;
            % 级联型
            % y[n] = a0*w[n] + a1*w[n-1] + a2*w[n-2]
            % w[n] = x[n] - b1*w[n-1] - b2*w[n-2]
            elseif "kCanonical" == obj.parameters.biquadCalcType
                 wn = xn - obj.coeffArray(obj.b1) * obj.stateArray(obj.x_z1) - obj.coeffArray(obj.b2) * obj.stateArray(obj.x_z2);
                 yn = obj.coeffArray(obj.a0) * wn + obj.coeffArray(obj.a1) * obj.stateArray(obj.x_z1) + obj.coeffArray(obj.a2) * obj.stateArray(obj.x_z2);
                 % 状态更新
                 obj.stateArray(obj.x_z2) = obj.stateArray(obj.x_z1);
                 obj.stateArray(obj.x_z1) = wn;
            % 转置直接型
            elseif "kTransposeDirect" == obj.parameters.biquadCalcType
                wn = xn + obj.stateArray(obj.y_z1);
                yn = obj.coeffArray(obj.a0) * wn + obj.stateArray(obj.x_z1);
                % 状态更新
                obj.stateArray(obj.y_z1) = obj.stateArray(obj.y_z2) - obj.coeffArray(obj.b1) * wn;
                obj.stateArray(obj.y_z2) = -obj.coeffArray(obj.b2) * wn;
                obj.stateArray(obj.x_z1) = obj.stateArray(obj.x_z2) + obj.coeffArray(obj.a1) * wn;
                obj.stateArray(obj.x_z2) = obj.coeffArray(obj.a2) * wn;
            % 转置级联型
            elseif "kTransposeCanonical" == obj.parameters.biquadCalcType
                yn = obj.coeffArray(obj.a0) * xn + obj.stateArray(obj.x_z1);
                % 状态更新
                obj.stateArray(obj.x_z1) = obj.coeffArray(obj.a1) * xn - obj.coeffArray(obj.b1) * yn + obj.stateArray(obj.x_z2);
                obj.stateArray(obj.x_z2) = obj.coeffArray(obj.a2) * xn - obj.coeffArray(obj.b2) * yn;
            else
                yn = xn;
            end
        end
        
    end
end