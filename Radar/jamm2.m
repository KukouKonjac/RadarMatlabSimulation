%% 雷达抗干扰功能级建模：频带不匹配与通过系数仿真
clear; clc; close all;

% --- 1. 参数设置 ---
fs = 800e6;                 % 采样率
T = 50e-6;                  % 观测时长
t = 0:1/fs:T-1/fs;
L = length(t);

% 干扰信号参数
B_j = 20e6;                 % 干扰带宽 (40MHz)
f_j_center = 50e6;          % 干扰中心频率 (50MHz)
[b_j, a_j] = butter(4, [(f_j_center - B_j/2), (f_j_center + B_j/2)]/(fs/2), 'bandpass');
s_jammer = filter(b_j, a_j, randn(1, L) + 1i*randn(1, L));

% 雷达接收机参数
B_rx = 20e6;                % 雷达带宽 (20MHz)
f_radar_list = 50e6 : 1e6 : 110e6; % 扫描频率范围

% 预计算数据
jamming_efficiency = zeros(1, length(f_radar_list));
[pxx_jammer, f_axis] = periodogram(s_jammer, rectwin(L), L, fs, 'centered');
pxx_jammer_norm = pxx_jammer / max(pxx_jammer); % 线性功率归一化

% --- 2. 仿真计算循环 ---
for i = 1:length(f_radar_list)
    f_0 = f_radar_list(i);
    Wn = [f_0 - B_rx/2, f_0 + B_rx/2] / (fs/2);
    Wn = max(min(Wn, 0.99), 0.01); 
    [b_rx, a_rx] = butter(4, Wn, 'bandpass');
    
    % 计算通过滤波器的干扰功率
    s_received = filter(b_rx, a_rx, s_jammer);
    jamming_efficiency(i) = var(s_received);
end

% --- 3. 绘制图 1：频域窗口与干扰遮盖示意图 ---
figure('Color', 'w', 'Name', '频域遮盖关系');
hold on;

% 绘制更加明显的干扰信号区域（使用鲜艳的颜色和半透明填充）
fill_f = f_axis / 1e6;
fill_p = 10*log10(pxx_jammer_norm + 1e-10);
fill([fill_f, fliplr(fill_f)], [fill_p, -100*ones(size(fill_p))], [0.9 0.4 0.4], ...
    'FaceAlpha', 0.4, 'EdgeColor', [0.8 0 0], 'LineWidth', 1, 'DisplayName', '干扰信号能量域');

% 仅选择三个具有代表性的中心频率：50MHz(全中), 70MHz(边缘), 95MHz(脱离)
plot_fc = [50e6, 70e6, 95e6];
colors = [0 0.45 0.74; 0.85 0.33 0.1; 0.47 0.67 0.19]; % 选定三种颜色

for k = 1:length(plot_fc)
    Wn_plot = [plot_fc(k) - B_rx/2, plot_fc(k) + B_rx/2] / (fs/2);
    [b_p, a_p] = butter(4, Wn_plot, 'bandpass');
    [h, w] = freqz(b_p, a_p, 2048, fs);
    plot(w/1e6, 20*log10(abs(h) + 1e-10), 'Color', colors(k,:), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('雷达窗口 fc=%dMHz', plot_fc(k)/1e6));
end

title('雷达工作窗口与干扰频带的重叠关系 (dB尺度)');
xlabel('频率 (MHz)'); ylabel('幅度/增益 (dB)');
axis([0, 200, -60, 5]); 
grid on; legend('Location', 'northeast');

% --- 1. 计算新的归一化自变量 (考虑 Bj) ---
% 这里的 B_sum_half 代表两个频谱边缘相碰所需的距离
B_sum_half = (B_rx + B_j) / 2;
x_norm = abs(f_radar_list - f_j_center) / B_sum_half;

% --- 2. 线性功率因数 rho ---
rho = jamming_efficiency / max(jamming_efficiency);

% --- 3. 通用功能级拟合模型 ---
% 由于自变量已经考虑了带宽，xc 现在可以取一个接近 1 的常数
x_fine = linspace(0, max(x_norm), 200);
xc = 0.9;  % 经验值：在边缘接触附近功率下降一半
n = 5;     % 陡峭度，与滤波器阶数正相关
rho_fit = 1 ./ (1 + (x_fine ./ xc).^(2 * n));

% --- 4. 绘图 ---
figure('Color', 'w', 'Name', '通用带宽匹配模型');
plot(x_norm, rho, 'r.', 'MarkerSize', 15, 'DisplayName', '仿真数据 (含Bj影响)');
hold on;
plot(x_fine, rho_fit, 'b-', 'LineWidth', 2, 'DisplayName', '拟合模型');

grid on;
xlabel('归一化分离度 x = |\Delta f| / [(B_{rx}+B_j)/2]');
ylabel('干扰通过系数 \rho');
title('考虑干扰带宽的通用功能级模型');
legend;
