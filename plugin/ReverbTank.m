% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% ReverbTank类

classdef ReverbTank < handle
    properties(Constant, Access = private)
        NUM_BRANCHES = 4;
        NUM_CHANNELS = 2;
    end
    
    properties(Access = public)
        parameters % 可调参数结构体
        % --- density
        % --- apfDelayMax_mSec
        % --- apfDelayWeight_Pct
        % --- fixeDelayMax_mSec
        % --- fixeDelayWeight_Pct
        % --- preDelayTime_mSec
        % --- lpf_g
        % --- kRT
        % --- lowShelf_fc
        % --- lowShelfBoostCut_dB
        % --- highShelf_fc
        % --- highShelfBoostCut_dB
        % --- wetLevel_dB
        % --- dryLevel_dB
        preDelay % SimpleDelay类
        branchDelays % SimpleDelay类数组
        branchNestedAPFs
        branchLPFs
        shelvingFilters
        apfDelayWeight
        fixedDelayWeight
        sampleRate
    end
    
    methods
        % === 构造函数 ===
        function obj = ReverbTank()
            obj.parameters = struct();
            obj.parameters.density = "kThick";
            obj.parameters.apfDelayMax_mSec = 5.0;
            obj.parameters.apfDelayWeight_Pct = 100.0;
            obj.parameters.fixeDelayMax_mSec = 50.0;
            obj.parameters.fixeDelayWeight_Pct = 100.0;
            obj.parameters.preDelayTime_mSec = 0.0;
            obj.parameters.lpf_g = 0.0;
            obj.parameters.kRT = 0.0;
            obj.parameters.lowShelf_fc = 0.0;
            obj.parameters.lowShelfBoostCut_dB = 0.0;
            obj.parameters.highShelf_fc = 0.0;
            obj.parameters.highShelfBoostCut_dB = 0.0;
            obj.parameters.wetLevel_dB = -3.0;
            obj.parameters.dryLevel_dB = -3.0;

            obj.preDelay = SimpleDelay();    
            
            obj.branchDelays = SimpleDelay.empty();
            obj.branchNestedAPFs = NestedDelayAPF.empty();
            obj.branchLPFs = SimpleLPF.empty();
            
            for i = 1 : obj.NUM_BRANCHES
                obj.branchDelays(i) = SimpleDelay();
                obj.branchNestedAPFs(i) = NestedDelayAPF();
                obj.branchLPFs(i) = SimpleLPF();
            end
            
            obj.shelvingFilters = TwoBandShelvingFilter.empty();
            for i = 1 : obj.NUM_CHANNELS
                obj.shelvingFilters(i) = TwoBandShelvingFilter();
            end
            
            obj.apfDelayWeight = [0.317, 0.873, 0.477, 0.291, 0.993, 0.757, 0.179, 0.575];
            obj.fixedDelayWeight = [1.0, 0.873, 0.707, 0.667];
            
            obj.sampleRate = 0.0;
        end
        
        % === 重置函数 ===
        function reset(obj, sampleRate_)
            obj.sampleRate = sampleRate_;
            
            obj.preDelay.reset(sampleRate_);
            obj.preDelay.createDelayBuffer(sampleRate_, 100.0);
            
            for i = 1 : obj.NUM_BRANCHES
                obj.branchDelays(i).reset(sampleRate_);
                obj.branchDelays(i).createDelayBuffer(sampleRate_, 100.0);
                
                obj.branchNestedAPFs(i).reset(sampleRate_);
                obj.branchNestedAPFs(i).createDelayBuffers(sampleRate_, 100.0, 100.0);
                
                obj.branchLPFs(i).reset(sampleRate_);
            end
            
            for i = 1 : obj.NUM_CHANNELS
                obj.shelvingFilters(i).reset(sampleRate_);
            end
            
        end
        
        function yn = processAudioSample(obj, xn)
            inputFrame = [xn, 0];
            outputFrame = obj.processAudioFrame(inputFrame, 1, 1);
            yn = outputFrame(1);
        end
        
        % === 数据处理函数 ===
        function outputFrame = processAudioFrame(obj, inputFrame, inputChannels, outputChannels)
            globFB = obj.branchDelays(obj.NUM_BRANCHES).readDelay();
            fb = obj.parameters.kRT * globFB;
            
            xnL = inputFrame(1);
            if inputChannels > 1
                xnR = inputChannels(2);
            else
                xnR = 0.0;
            end
            monoXn = 1.0 / inputChannels * xnL + 1.0 / inputChannels * xnR;
            
            preDelayOut = obj.preDelay.processAudioSample(monoXn);
            
            input = preDelayOut + fb;
            for i = 1 : obj.NUM_BRANCHES
                apfOut = obj.branchNestedAPFs(i).processAudioSample(input);
                lpfOut = obj.branchLPFs(i).processAudioSample(apfOut);
                delayOut = obj.parameters.kRT * obj.branchDelays(i).processAudioSample(lpfOut);
                input = delayOut + preDelayOut;
            end
            
            weight = 0.707;
            
            outL = 0.0;
            outL = outL + weight * obj.branchDelays(1).readDelayAtPercentage(23.0);
            outL = outL - weight * obj.branchDelays(2).readDelayAtPercentage(41.0);
            outL = outL + weight * obj.branchDelays(3).readDelayAtPercentage(59.0);
            outL = outL - weight * obj.branchDelays(4).readDelayAtPercentage(73.0);
            
            outR = 0.0;
            outR = outR - weight * obj.branchDelays(1).readDelayAtPercentage(29.0);
            outR = outR + weight * obj.branchDelays(2).readDelayAtPercentage(43.0);
            outR = outR - weight * obj.branchDelays(3).readDelayAtPercentage(61.0);
            outR = outR + weight * obj.branchDelays(4).readDelayAtPercentage(79.0);
            
            if "kThick" == obj.parameters.density
                outL = outL + weight*obj.branchDelays(1).readDelayAtPercentage(31.0);
                outL = outL - weight*obj.branchDelays(2).readDelayAtPercentage(47.0);
                outL = outL + weight*obj.branchDelays(3).readDelayAtPercentage(67.0);
                outL = outL - weight*obj.branchDelays(4).readDelayAtPercentage(83.0);

                outR = outR - weight*obj.branchDelays(1).readDelayAtPercentage(37.0);
                outR = outR + weight*obj.branchDelays(2).readDelayAtPercentage(53.0);
                outR = outR - weight*obj.branchDelays(3).readDelayAtPercentage(71.0);
                outR = outR + weight*obj.branchDelays(4).readDelayAtPercentage(89.0);
            end
            
            tankOutL = obj.shelvingFilters(1).processAudioSample(outL);
            tankOutR = obj.shelvingFilters(2).processAudioSample(outR);
            
            dry = 10 ^ (obj.parameters.dryLevel_dB / 20);
            wet = 10 ^ (obj.parameters.wetLevel_dB / 20);
            
            if 1 == outputChannels
                outputFrame(1) = dry * xnL + wet * (0.5 * tankOutL + 0.5 * tankOutR);
            else
                outputFrame(1) = dry * xnL + wet * tankOutL;
                outputFrame(2) = dry * xnR + wet * tankOutR;
            end
        end
        
        % === 获取参数 ===
        function parameters = getParameters(obj)
            parameters = obj.parameters;
        end
        
        % === 设置参数 ===
        function setParameters(obj, params)
            filterParams = obj.shelvingFilters(1).getParameters();
            filterParams.highShelf_fc = params.highShelf_fc;
            filterParams.highShelfBoostCut_dB = params.highShelfBoostCut_dB;
            filterParams.lowShelf_fc = params.lowShelf_fc;
            filterParams.lowShelfBoostCut_dB = params.lowShelfBoostCut_dB;
            
            obj.shelvingFilters(1).setParameters(filterParams);
            obj.shelvingFilters(2).setParameters(filterParams);
            
            lpfParams = obj.branchLPFs(1).getParameters();
            lpfParams.g = params.lpf_g;
            
            for i = 1 : obj.NUM_BRANCHES
                obj.branchLPFs(i).setParameters(lpfParams);
            end
            
            % 更新延迟模块
            delayParams = obj.preDelay.getParameters();
            delayParams.delayTime_mSec = params.preDelayTime_mSec;
            obj.preDelay.setParameters(delayParams);
            
            % 更新NestedAPF模块
            apfParams = obj.branchNestedAPFs(1).getParameters();
            delayParams = obj.branchDelays(1).getParameters();
            
            globalAPFMaxDelay = (obj.parameters.apfDelayWeight_Pct / 100) * obj.parameters.apfDelayMax_mSec;
            globalFixedMaxDelay = (obj.parameters.fixeDelayWeight_Pct / 100) * obj.parameters.fixeDelayMax_mSec;
            
            apfParams.enableLFO = true;
            apfParams.lfoMaxModulation_mSec = 0.3;
            apfParams.lfoDepth = 1.0;

            m = 1;
            for i = 1 : obj.NUM_BRANCHES
                apfParams.outerAPFdelayTime_mSec = globalAPFMaxDelay * obj.apfDelayWeight(m);
                m = m + 1;
                apfParams.innerAPFdelayTime_mSec = globalAPFMaxDelay * obj.apfDelayWeight(m);
                m = m + 1;
                apfParams.innerAPF_g = -0.5;
                apfParams.outerAPF_g = 0.5;
                
                if 1 == i
                    apfParams.lfoRate_Hz = 0.15;
                elseif 2 == i
                    apfParams.lfoRate_Hz = 0.33;
                elseif 3 == i
                    apfParams.lfoRate_Hz = 0.57;
                elseif 4 == i
                    apfParams.lfoRate_Hz = 0.73;
                end
                
                obj.branchNestedAPFs(i).setParameters(apfParams);
                
                delayParams.delayTime_mSec = globalFixedMaxDelay * obj.fixedDelayWeight(i);
                obj.branchDelays(i).setParameters(delayParams);
            end
            
            obj.parameters = params;
            
        end    
        
    end
end
