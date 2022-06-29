%% 一阶二阶滤波器仿真
% 参考<Audio Effects theory,implementation and application>
clc; clear all; close all;

%% N阶低通滤波器和高通滤波器
w = 0 : pi/50: pi;
z = exp(1i*w);
Nmax = 2; % 仿真最大阶数

% 标准滤波器
figure()
for N = 1:Nmax
    H = ones(1, length(w));
    for n = 1:N
        gamma = ((2*n - 1)/N - 1) * pi / 4;
        k = 1 / (2 * cos(gamma));
        q = -1;
        p = 1i * tan(gamma);
        Hall = k * (z - q)./(z - p); % (3.5)
        H = H.*Hall;
    end
    subplot(2,2,n);
    plot(w/pi, abs(H));
    title(sprintf("%d阶低通滤波器", n));
end

% 改变截至频率位置
% 利用映射函数F = (z - beta) ./ (1 - beta*z)
wc = pi/8;
beta = (1 - tan(wc/2)) / (1 + tan(wc/2)); % (3.15)
znew = (z - beta)./(1 - beta*z); % (3.14)

figure()
for N = 1:Nmax
    H = ones(1, length(w));
    for n = 1:N
        gamma = ((2*n - 1)/N - 1) * pi / 4;
        k = 1 / (2 * cos(gamma));
        q = -1;
        p = 1i * tan(gamma);
        Hall = k * (znew - q)./(znew - p); % (3.16)
        H = H.*Hall;
    end
    subplot(2,2,n);
    plot(w/pi, abs(H));
    title(sprintf("%d阶低通滤波器", n));
end

% 改为高通滤波器
% 利用映射函数F = -z，这个映射函数实际上是把相位平移了pi，这一点可以通过单位圆进行理解
zhp = -z;

figure()
for N = 1:Nmax
    H = ones(1, length(w));
    for n = 1:N
        gamma = ((2*n - 1)/N - 1) * pi / 4;
        k = 1 / (2 * cos(gamma));
        q = -1;
        p = 1i * tan(gamma);
        Hall = k * (zhp - q)./(zhp - p); % (3.16)
        H = H.*Hall;
    end
    subplot(2,2,n);
    plot(w/pi, abs(H));
    title(sprintf("%d阶高通滤波器", n));
end

% 高通滤波器截至频率变换
% 先变换低通滤波器的截至评频率，再平移pi
wc = pi/8;

beta = (1 - tan(wc/2)) / (1 + tan(wc/2)); % (3.15)
znew = (z - beta)./(1 - beta*z); % (3.14)
zhpnew = - znew;

figure()
for N = 1:Nmax
    H = ones(1, length(w));
    for n = 1:N
        gamma = ((2*n - 1)/N - 1) * pi / 4;
        k = 1 / (2 * cos(gamma));
        q = -1;
        p = 1i * tan(gamma);
        Hall = k * (zhpnew - q)./(zhpnew - p); % (3.16)
        H = H.*Hall;
    end
    subplot(2,2,n);
    plot(w/pi, abs(H));
end  

%% N阶low-shelf滤波器
w = 0 : pi/50: pi;
z = exp(1i*w);
Nmax = 4; % 仿真最大阶数
gmax = 4; %最大增益

% 标准滤波器
figure()
for N = 1:Nmax
    H = ones(1, length(w));
    for n = 1:N
        gamma = ((2*n - 1)/N - 1) * pi / 4;
        k = 1 / 2;
        g = gmax^(1/N);
        p = 1i * tan(gamma);
        Hall = k *...
               ((1 + p + g*(1 - p))*z - (1 + p - g*(1 - p)))./...
               (z - p); % (3.18)
        H = H.*Hall;
    end
    subplot(2,2,n);
    plot(w/pi, abs(H));
    title(sprintf("%d阶low-shelf滤波器", n));
    axis([0, 1, 0, gmax]);
end

% 改变截至频率位置
% 利用映射函数F = (z - beta) / (1 - beta*z)
wc = pi/8;
beta = (1 - tan(wc/2)) / (1 + tan(wc/2)); % (3.15)
znew = (z - beta)./(1 - beta*z); % (3.14)

figure()
for N = 1:Nmax
    H = ones(1, length(w));
    for n = 1:N
        gamma = ((2*n - 1)/N - 1) * pi / 4;
        k = 1 / 2;
        g = gmax^(1/N);
        p = 1i * tan(gamma);
        Hall = k *...
               ((1 + p + g*(1 - p))*znew - (1 + p - g*(1 - p)))./...
               (znew - p); % (3.18)
        H = H.*Hall;
    end
    subplot(2,2,n);
    plot(w/pi, abs(H));
    title(sprintf("%d阶low-shelf滤波器", n));
    axis([0, 1, 0, gmax]);
end

%% 带通滤波器
% 首先构建低通滤波器原型
% 然后利用映射函数 F = -(z.^2 + a1 * z + a2) ./ (a2*z.^2 + a1*)
w = 0 : pi/50: pi;
z = exp(1i*w);
Nmax = 4; % 仿真最大阶数
wc = pi/4; % 带通滤波器中心频率
B = pi/8;

% 标准滤波器
figure()
for N = 1:Nmax
    H = ones(1, length(w));
    for n = 1:N
        gamma = ((2*n - 1)/N - 1) * pi / 4;
        k = 1 / (2 * cos(gamma));
        q = -1;
        p = 1i * tan(gamma);
        Hall = k *...
               (z.^2 - z*(1 + q)*cos(wc) + q)./ ...
               (z.^2 - z*(1 + p)*cos(wc) + p); % (3.21)
        H = H.*Hall;
    end
    subplot(2,2,n);
    plot(w/pi, abs(H));
    title(sprintf("%d阶低通滤波器", n));
end


% 考虑带宽
% 先将低通截止频率移动到带宽处
% 再将其改为带通
% 注意低通滤波器截止频率改变后，(3.21)式子用到的p和q都要改变
beta = (1 - tan(B/2)) / (1 + tan(B/2)); % (3.15)
figure()
for N = 1:Nmax
    H = ones(1, length(w));
    for n = 1:N
        gamma = ((2*n - 1)/N - 1) * pi / 4;
        k = 1 / (2 * cos(gamma));
        q = -1;
        p = 1i * tan(gamma);
        qnew = (q + beta)/(1 + beta*q); % 映射得到新值
        pnew = (p+beta)/(1 + beta*p); % 映射得到新知
        Hall = k *...
               (z.^2 - z*(1 + qnew)*cos(wc) + qnew)./ ...
               (z.^2 - z*(1 + pnew)*cos(wc) + pnew); 
        H = H.*Hall;
    end
    subplot(2,2,n);
    plot(w/pi, abs(H));
    title(sprintf("%d阶低通滤波器", n));
end

% 考虑带宽
% 先将低通截止频率移动到带宽处
% 再将其改为带通
% 注意低通滤波器截止频率改变后，(3.21)式子用到的p和q都要改变
beta = (1 - tan(B/2)) / (1 + tan(B/2)); % (3.15)
znew = (z - beta)./(1 - beta*z); % (3.14)

figure()
for N = 1:Nmax
    H = ones(1, length(w));
    for n = 1:N
        gamma = ((2*n - 1)/N - 1) * pi / 4;
        k = 1 / (2 * cos(gamma));
        q = -1;
        p = 1i * tan(gamma);
        qnew = (q + beta)/(1 + beta*q); % 映射得到新值
        pnew = (p+beta)/(1 + beta*p); % 映射得到新知
        Hall = k *...
               (znew.^2 - znew*(1 + q)*cos(wc) + q)./ ...
               (znew.^2 - znew*(1 + p)*cos(wc) + p); 
        H = H.*Hall;
    end
    subplot(2,2,n);
    plot(w/pi, abs(H));
    title(sprintf("%d阶低通滤波器", n));
end

%% table 3.1
%% 一阶低通
wc = pi/4;
b = [tan(wc/2), tan(wc/2)];
a = [1+tan(wc/2), tan(wc/2)-1];
[h, w] = freqz(b, a, 'whole', 2001);
figure;
subplot(2,1,1);
plot(w/pi,abs(h),'.-')
axis([0 1 0 1])
legend('Response','Passband','Location','South')
ylabel('Magnitude')
title("一阶低通")
subplot(2,1,2);
plot(w/pi,db(h),'.-')
axis([0 1 -60 10])
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')

%% 一阶高通
wc = pi/4;
b = [1, -1];
a = [1+tan(wc/2), tan(wc/2)-1];
[h, w] = freqz(b, a, 'whole', 2001);
figure;
plot(w/pi,abs(h),'.-')
axis([0 1 -1 2])
legend('Response','Passband','Location','South')
ylabel('Magnitude')
title("一阶高通")

%% 一阶low-shelf
wc = pi/4;
G = 2;
b = [1+G*tan(wc/2), G*tan(wc/2)-1];
a = [1+tan(wc/2), tan(wc/2)-1];
[h, w] = freqz(b, a, 'whole', 2001);
figure;
plot(w/pi,abs(h),'.-')
axis([0 1 -1 2])
legend('Response','Passband','Location','South')
ylabel('Magnitude')
title("一阶low-shelf")

%% 一阶high-shelf
wc = pi/4;
G = 2;
b = [tan(wc/2)+G, tan(wc/2)-G];
a = [1+tan(wc/2), tan(wc/2)-1];
[h, w] = freqz(b, a, 'whole', 2001);
figure;
plot(w/pi,abs(h),'.-')
axis([0 1 -1 2])
legend('Response','Passband','Location','South')
ylabel('Magnitude')
title("一阶high-shelf")

%% 二阶带通
wc = pi/4;
B = pi/8;
b = [tan(B/2), 0, -tan(B/2)];
a = [1+tan(B/2), -2*cos(wc), 1-tan(B/2)];
[h, w] = freqz(b, a, 'whole', 2001);
figure;
plot(w/pi,abs(h),'.-')
axis([0 1 -1 2])
legend('Response','Passband','Location','South')
ylabel('Magnitude')
title("二阶带通")

%% 二阶带阻
wc = pi/4;
B = pi/8;
b = [1, -2*cos(wc), 1];
a = [1+tan(B/2), -2*cos(wc), 1-tan(B/2)];
[h, w] = freqz(b, a, 'whole', 2001);
figure;
plot(w/pi,abs(h),'.-')
axis([0 1 -1 2])
legend('Response','Passband','Location','South')
ylabel('Magnitude')
title("二阶带阻")

%% Peaking/Notch
wc = pi/4;
G = 2;
B = pi/8;
b = [1+G*tan(B/2), -2*cos(wc), 1-G*tan(B/2)];
a = [1+tan(B/2), -2*cos(wc), 1-tan(B/2)];
[h, w] = freqz(b, a, 'whole', 2001);
figure;
plot(w/pi,abs(h),'.-')
axis([0 1 -1 2])
legend('Response','Passband','Location','South')
ylabel('Magnitude')
title("Peaking/Notch")





