%%
clc; clear all; close all;

%%
[audio, fs] = audioread("data/0.wav");
gain=0.8;
simpledelay_output = zeros(length(audio), 1);

b = [1 1 1 1];
c = [0.8 0.8 0.8 0.8];
% Feedback matrix
% 反馈矩阵是一个酉矩阵
a(1, : ) =[0 1 1 0];
a(2, : ) =[-1 0 0 -1];
a(3, : ) =[1 0 0 -1];
a(4, : ) =[0 1 -1 0];
a2=a*(1/sqrt(2) ) * gain;

% Delay lines, use prime numbers
% 使用素数，防止不同延时之间产生相干性，导致结果染色
m=[149 211 263 293]';
z1=zeros(1, max(max(m) ) );
z2=zeros(1, max(max(m) ) );
z3=zeros(1, max(max(m) ) );
z4=zeros(1, max(max(m) ) );


for n = 1: length(audio)
    tmp = [z1(m(1) ) z2(m(2) ) z3(m(3) ) z4(m(4) ) ];
    simpledelay_output(n) = audio(n) + c(1) *z1(m(1) ) + c(2) *z2(m(2) ) ...
    + c(3) *z3(m(3) ) + c(4) *z4(m(4) );
    z1 = [(audio(n)*b(1) + tmp*a2(1, : )') z1(1: length(z1) -1)];
    z2 = [(audio(n)*b(2) + tmp*a2(2, : )') z2(1: length(z2) -1)];
    z3 = [(audio(n)*b(3) + tmp*a2(3, : )') z3(1: length(z3) -1)];
    z4 = [(audio(n)*b(4) + tmp*a2(4, : )') z4(1: length(z4) -1)];
end

audiowrite('data/FDN_0.wav', simpledelay_output, fs);