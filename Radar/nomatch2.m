% 定义参数
N = 128; % 滤波器长度
fs = 1000; % 采样频率
f0 = 50; % 通带中心频率
f1 = 150; % 阻带中心频率
BW = (f1 - f0) / fs; % 带通宽度
% 生成汉宁窗
h = hamming(N);
% 设计滤波器
[b, a] = designfilt('fir1', N, f0, BW, fs);
% 生成带噪声信号
x = sin(2*pi*f0*t) + 0.5*sin(2*pi*f1*t) + randn(size(t));
y = filter(b, a, x);
% 绘图
figure;
subplot(2, 1, 1); plot(t, x); title('Input Signal'); xlabel('Time'); ylabel('Amplitude');
subplot(2, 1, 2); plot(t, y); title('Filtered Signal'); xlabel('Time'); ylabel('Amplitude');aw