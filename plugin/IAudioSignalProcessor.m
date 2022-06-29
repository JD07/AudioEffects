% 参考Designing Audio Effect Plugins in C++的fxobjects文件，实现数字滤波器的各个基本组件
% 基础类，给子类提供函数模板

classdef IAudioSignalProcessor < handle
    methods(Access = public)
        function obj = IAudioSignalProcessor()
        end
                
        function reset(obj, sampleRate_)
        end

        function processAudioSample(obj, ~)
        end

        function flag = canProcessAudioFrame(obj)
            flag = false;
            return;
        end
        
        function setSampleRate(obj, sampleRate_)
        end
        
        function enableAuxInput(obj, enableAuxInput)
        end
        

        function yn = processAuxInputAudioSample(xn)
            yn = xn;
            return;
        end

        function outputFrame = processAudioFrame(obj, inputFrame, inputChannels, outputChannels)
            outputFrame = inputFrame;
            return;
        end
    end
end