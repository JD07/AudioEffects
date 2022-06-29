% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% 低频振荡器模块
% 使用方法：
% （1）调用构造函数初始化类的实例
% （2）调用reset方法设置采样率
% （3）调用getParameters方法，返回可调整参数结构体
% （4）设置可调整参数结构体，并利用setParameters方法设置参数

% 2021/2/21 通过ModuleTest测试

classdef LFO < handle
    properties(Access = public)
        lfoParameters % 可调整参数结构体
        % --- waveform % 波形: kSin kTriangle kSaw
        % --- frequency_Hz % 振荡器频率
        sampleRate % 采样率
        modCounter % modulo counter [0.0, +1.0]
        phaseInc % 相位步长 phaseInc = fo/fs
        modCounterQP % Quad Phase modulo counter [0.0, +1.0]
    end
    
    methods(Access = public)
        % === 构造函数 ===
        % 按照默认参数初始化
        function obj = LFO()
            obj.lfoParameters = struct();
            obj.lfoParameters.waveform = "kSin";
            obj.lfoParameters.frequency_Hz = 0;

            obj.sampleRate = 0.0;
            obj.modCounter = 0.0;
            obj.phaseInc = 0.0;
            obj.modCounterQP = 0.25;
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            obj.sampleRate = sampleRate_;
            obj.phaseInc = obj.lfoParameters.frequency_Hz / obj.sampleRate;
            
            %  时间基准变量
            obj.modCounter = 0.0; % < modulo counter [0.0, +1.0]
            obj.modCounterQP = 0.25; % < Quad Phase modulo counter [0.0, +1.0]
        end
        
        % === 参数结构体返回函数 ===
        function lfoParameters = getParameters(obj)
            lfoParameters = obj.lfoParameters;
        end
        
        % === 参数设置函数 ===
        function setParameters(obj, params)
            obj.lfoParameters = params;
            obj.phaseInc = params.frequency_Hz / obj.sampleRate;
        end
          
        % === 渲染新的音频输出结构体 ===
        function output = renderAudioOutput(obj)
            % 音频输出结构体初始化
            output = struct();
            output.normalOutput = 0.0;
            output.invertedOutput = 0.0;
            output.quadPhaseOutput_pos = 0.0;
            output.quadPhaseOutput_neg = 0.0;
            
            % 首先检查计数器
            obj.modCounter = checkAndWrapModulo(obj.modCounter, obj.phaseInc);
            % QP输出总时跟随当前模计数器，先令两者相等
            obj.modCounterQP = obj.modCounter;
            % 然后将QP计数器前进0.25相位即90°
            obj.modCounterQP = advanceAndCheckWrapModulo(obj.modCounterQP, 0.25);
            
            % 计算振荡器数值
            waveform = obj.lfoParameters.waveform;
            if waveform == "kSin"
                % calculate normal angle
                % 我们使用抛物线函数模拟正弦曲线，而该函数的输入范围为[-pi pi]
                % 所以我们需要将角度范围转化为双极性
                angle = obj.modCounter*2.0*pi - pi;
                % norm output with obj.parabolicSine approximation
                output.normalOutput = parabolicSine(-angle);
                % calculate QP angle
                angle = obj.modCounterQP*2.0*pi - pi;
                % calc QP output
                output.quadPhaseOutput_pos = parabolicSine(-angle);
            elseif waveform == "kTriangle"
                % triv saw
                output.normalOutput = unipolarToBipolar(obj.modCounter);
                % bipolar triagle
                output.normalOutput = 2.0 * abs(output.normalOutput) - 1.0;
                % quad phase
                output.quadPhaseOutput_pos = unipolarToBipolar(obj.modCounterQP);
                % bipolar triagle
                output.quadPhaseOutput_pos = 2.0 * abs(output.quadPhaseOutput_pos) - 1.0;
            elseif waveform == "kSaw"
                output.normalOutput = unipolarToBipolar(obj.modCounter);
                output.quadPhaseOutput_pos = unipolarToBipolar(obj.modCounterQP);
            end
            output.quadPhaseOutput_neg = - output.quadPhaseOutput_pos;
            output.invertedOutput = - output.normalOutput;
            
            obj.modCounter = advanceModulo(obj.modCounter, obj.phaseInc);
        end

        % === 根据输入的偏移，计算波形输出
        function yn = plusOffset(obj, offset)
            % 首先检查计数器
            obj.modCounter = checkAndWrapModulo(obj.modCounter, obj.phaseInc);
            % 先令偏移计数器和lfo计数器相等
            modCounterOffset = obj.modCounter;
            % 然后令偏移计数器加上偏移
            modCounterOffset = advanceAndCheckWrapModulo(modCounterOffset, offset);
            
            % 计算振荡器数值
            waveform = obj.lfoParameters.waveform;
            if waveform == "kSin"
                % calculate normal angle
                % 我们使用抛物线函数模拟正弦曲线，而该函数的输入范围为[-pi pi]
                % 所以我们需要将角度范围转化为双极性
                angle = modCounterOffset*2.0*pi - pi;
                % norm output with obj.parabolicSine approximation
                yn = parabolicSine(-angle);
            elseif waveform == "kTriangle"
                % triv saw
                yn = unipolarToBipolar(modCounterOffset);
                % bipolar triagle
                yn = 2.0 * abs(yn) - 1.0;
            elseif waveform == "kSaw"
                yn = unipolarToBipolar(modCounterOffset);
            end
        end
    end
    
end

% === 检查模计数器，并在需要的时候进行调整 ===
function moduloCounter = checkAndWrapModulo(moduloCounter, phaseInc)
    % 对于正频率
    if phaseInc > 0 && moduloCounter >= 1.0
        moduloCounter = moduloCounter - 1.0;
    end
    % 对于负频率
    if phaseInc < 0 && moduloCounter <= 0.0
        moduloCounter = moduloCounter + 1.0;
    end
end

% === 模计数器计数 ===
function moduloCounter = advanceModulo(moduloCounter, phaseInc)
    moduloCounter = moduloCounter + phaseInc;
end

% === 模计数器计数,然后检查计数器，并在需要的时候进行调整 ===
function moduloCounter = advanceAndCheckWrapModulo(moduloCounter, phaseInc)
    moduloCounter = moduloCounter + phaseInc;
    % 对于正频率
    if phaseInc > 0 && moduloCounter >= 1.0
        moduloCounter = moduloCounter - 1.0;
    end
    % 对于负频率
    if phaseInc < 0 && moduloCounter <= 0.0
        moduloCounter = moduloCounter + 1.0;
    end
end

% === 将单极性转为双极性 ===
function output =  unipolarToBipolar(value)
    output = 2.0*value - 1.0;
end

% === 利用抛物线近似正弦曲线 ===
% 注意，输入范围为-pi~pi
function y = parabolicSine(angle)
    B = 1.2732; % 4.0 / pi
    C = -0.4053; % -4.0 / pi
    P = 0.225;

    y = B * angle + C * angle * abs(angle);
    y = P * (y * abs(y) - y) + y;
end
