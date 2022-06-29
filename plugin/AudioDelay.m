% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% 音频延时模块
% 三种处理模式：单输入单输出，单输入双输出，双输入双输出
% delayAlgorithm { kNormal, kPingPong }
% delayUpdateType { kLeftAndRight, kLeftPlusRatio } ...
% kLeftAndRoght使用输入的左右延时值来初始化，而kLeftPlusRatio则是左声道使用输入值初始化...
% 右声道则根据这个值进行计算，从而确保两边声道的延时成特定的比值
% Ping-Pong需要两个通道的延时成比例才能起到效果。
% 使用方法：
% （1）调用构造函数初始化类的实例
% （2）调用reset方法设置采样率
% （3）调用getParameters方法，返回可调整参数结构体
% （4）设置可调整参数结构体，并利用setParameters方法设置参数

% 2021/2/22 测试未发现异常
% 2021/2/23 改为handle类，测试未发现异常

classdef AudioDelay < handle
    properties(Access = public)
        parameters % 可调整参数结构体
        % --- algorithm % 使用的算法
        % --- wetLevel_dB % 湿信号等级
        % --- dryLevel_dB % 干信号等级
        % --- feedback_Pct % 反馈比例
        % --- updateType % 左右声道延时的更新方式
        % --- leftDelay_mSec % 左声道延时
        % --- rightDelay_mSec % 右声道延时
        % --- delayRatio_Pct % 延时比例，范围0-100
        % --- reverseChannels % 是否将输出相互交叉
        % --- interpolate % 插值方法： kNearest, kLinear, kSquare, kCubic
        
        sampleRate % 采样率
        samplesPerMSec % 每ms对应的采样点个数
        delayInSamples_L % 左声道延时点数
        delayInSamples_R % 右声道延时点数
        bufferLength_mSec % 缓存区时间长度
        bufferLength % 缓存区对应的点数
        wetMix % 湿信号系数
        dryMix % 干信号系数

        delayBuffer_L % 类变量，左声道缓存区
        delayBuffer_R % 类变量，右声道缓存区
    end
    
    methods(Access = public)
        % === 构造函数 ===
        % 按照默认参数初始化
        function obj = AudioDelay()
            % 可调参数结构体初始化
            obj.parameters = struct();
            obj.parameters.algorithm = "kNormal";
            obj.parameters.wetLevel_dB = -3.0;
            obj.parameters.dryLevel_dB = -3.0;
            obj.parameters.feedback_Pct = 0.0;
            obj.parameters.updateType = "kLeftAndRight";
            obj.parameters.leftDelay_mSec = 0.0;
            obj.parameters.rightDelay_mSec = 0.0;
            obj.parameters.delayRatio_Pct = 100.0;
            obj.parameters.reverseChannels = false;
            obj.parameters.interpolate = "kNearest";
            
            % 参数初始化
            obj.sampleRate = 0.0;
            obj.samplesPerMSec = 0.0;
            obj.delayInSamples_L = 0.0;
            obj.delayInSamples_R = 0.0;
            obj.bufferLength_mSec = 0.0;
            obj.bufferLength = 0;
            obj.wetMix = 0.707; % 默认-3dB
            obj.dryMix = 0.707; % 默认-3dB
            obj.delayBuffer_L = CircularBuffer();
            obj.delayBuffer_L.setInterpolate(obj.parameters.interpolate);
            obj.delayBuffer_R = CircularBuffer();
            obj.delayBuffer_R.setInterpolate(obj.parameters.interpolate);
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate)
            if sampleRate == obj.sampleRate 
                % 如果采样率没有发生变化，则只需要清空延迟缓存区
                obj.delayBuffer_L.flushBuffer();
                obj.delayBuffer_R.flushBuffer();
            else
                % 采样率发生变化，需要保存新的采样率并创建新的延时缓存
                obj.createDelayBuffers(sampleRate, obj.bufferLength_mSec);
            end
        end
        
        % === 获取参数 ===
        function parameters = getParameters(obj)
            parameters = obj.parameters;
        end

        % === 参数设置函数 ===
        function setParameters(obj, params)
            if params.dryLevel_dB ~= obj.parameters.dryLevel_dB
                obj.dryMix = 10^(params.dryLevel_dB/20);
            end
            if params.wetLevel_dB ~= obj.parameters.wetLevel_dB
                obj.wetMix = 10^(params.wetLevel_dB/20);
            end
            obj.parameters = params;
            
            % 若延时时间大于缓存区长度，则重新分配内存
            maxDelay_mSec = max(obj.parameters.leftDelay_mSec, obj.parameters.rightDelay_mSec);
            if maxDelay_mSec > obj.bufferLength_mSec
                obj.createDelayBuffers(obj.sampleRate, maxDelay_mSec);
            end
            
            % 根据配置设置剩余参数
            if obj.parameters.updateType == "kLeftAndRight"
                % 根据输入的左右声道延时值设置左右声道
                newDelayInSamples_L = obj.parameters.leftDelay_mSec * obj.samplesPerMSec;
                newDelayInSamples_R = obj.parameters.rightDelay_mSec * obj.samplesPerMSec;
                obj.delayInSamples_L = newDelayInSamples_L;
                obj.delayInSamples_R = newDelayInSamples_R;
            elseif obj.parameters.updateType == "kLeftPlusRatio"
                % 根据输入的左声道延时值设置左声道，然后更具delayRatio设置右声道
                % Ping-Pong延时要求左右声道的延时成比例，不能一样
                delayRatio = obj.parameters.delayRatio_Pct / 100.0;
                if delayRatio < 0.0
                    delayRatio = 0.0;
                elseif delayRatio > 1.0
                    delayRatio = 1.0;
                end
                
                newDelayInSamples = obj.parameters.leftDelay_mSec * obj.samplesPerMSec;
                
                obj.delayInSamples_L = newDelayInSamples;
                obj.delayInSamples_R = obj.delayInSamples_L * delayRatio;
            end
            obj.delayBuffer_L.setInterpolate(obj.parameters.interpolate);
            obj.delayBuffer_R.setInterpolate(obj.parameters.interpolate);
        end
        
        % 延时缓存创造函数
        function createDelayBuffers(obj, sampleRate_, bufferLength_mSec_)
            obj.bufferLength_mSec = bufferLength_mSec_;
            obj.sampleRate = sampleRate_;
            obj.samplesPerMSec = obj.sampleRate / 1000.0;

            obj.bufferLength = ceil(obj.bufferLength_mSec * obj.samplesPerMSec);
            obj.delayBuffer_L.createCircularBuffer(obj.bufferLength);
            obj.delayBuffer_R.createCircularBuffer(obj.bufferLength);
        end
        
        % 处理单通道数据
        function output = processAudioSample(obj, xn)
            yn = obj.delayBuffer_L.readBuffer(obj.delayInSamples_L);
            dn = xn + (obj.parameters.feedback_Pct/100) * yn;
            obj.delayBuffer_L.writeBuffer(dn);
            output = obj.dryMix * xn + obj.wetMix * yn;
        end
        
        % 处理单帧数据
        function outputFrame = processAudioFrame(obj, inputFrame, inputChannels, outputChannels)
            if inputChannels == 0 || outputChannels == 0
                error('Error : InputChannels and outputChannels should not be 0');
            end
            if "kNormal" ~= obj.parameters.algorithm && "kPingPong" ~= obj.parameters.algorithm
                error('Error : Algorithm should be kNormal or kPingPong');
            end
            % 如果只输出单通道，转到单输入单输出
            if outputChannels == 1
                output = obj.processAudioSample(inputFrame(1));
                outputFrame(1) = output;
                return;
            end
            
            % 如果输出双通道
            % 左通道
            xnL = inputFrame(1);
            % 右通道（如果是单声道输入，则复制左声道的值）
            if inputChannels > 1
                xnR = inputFrame(2);
            else
                xnR = xnL;
            end
            % 读取左声道延迟值
            ynL = obj.delayBuffer_L.readBuffer(obj.delayInSamples_L);
            % 读取右声道延迟值
            ynR = obj.delayBuffer_R.readBuffer(obj.delayInSamples_R);
            % 计算左声道延迟单元输入值
            dnL = xnL + (obj.parameters.feedback_Pct / 100) * ynL;
            % 计算右声道延迟单元输入值
            dnR = xnR + (obj.parameters.feedback_Pct / 100) * ynR;
            
            % 解码
            if obj.parameters.algorithm == "kNormal"
                % 正常模式下，两侧对应延时
                obj.delayBuffer_L.writeBuffer(dnL);
                obj.delayBuffer_R.writeBuffer(dnR);
            elseif obj.parameters.algorithm == "kPingPong"
                % 乒乓模式下，两侧交叉延时
                obj.delayBuffer_L.writeBuffer(dnR);
                obj.delayBuffer_R.writeBuffer(dnL);
            end
            
            if obj.parameters.reverseChannels
                outputL = obj.dryMix * xnL + obj.wetMix * ynR;
                outputR = obj.dryMix * xnR + obj.wetMix * ynL;
            else
                outputL = obj.dryMix * xnL + obj.wetMix * ynL;
                outputR = obj.dryMix * xnR + obj.wetMix * ynR;
            end

            
            outputFrame(1) = outputL;
            outputFrame(2) = outputR;
        end
    end
end