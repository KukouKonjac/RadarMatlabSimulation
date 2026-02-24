close all; clear; clc;

%% ================= LFM 信号参数 =================
B = 2e6;
T = 100e-6;
u = B / T;
fs = 4 * B;
N = round(T * fs);
t = (0:N-1).' / fs;
B_rx = B;                 % 接收机等效带宽
f_cut = B_rx/2;           % 基带截止频率

bp_order = 256;
rx_lpf = fir1(bp_order, f_cut/(fs/2), 'low');



%% ================= 雷达与噪声参数 =================
Pt = 1e10;
Gt = 30; Gr = 30;
sigma = 1;
Lsys = 1;

k = physconst('Boltzmann');
T0 = 290;
NF = 3;
B_rx = fs;

c = physconst('LightSpeed');
fc_radar = 10e9;
lambda = c / fc_radar;

%% ================= 发射信号 =================
s_tx = exp(1j*pi*u*t.^2) ;

%% ================= 设置目标距离列表 =================
R_targets = [1500];  % 可以自己添加或修改多个距离
colors = lines(length(R_targets)); % 绘图颜色

%% ================= 设置频偏列表 =================
df_list = 0 : 0.05*B : 1.2*B;

%% ================= 绘制 SNR vs 频偏 =================
figure; hold on;
for r_idx = 1:length(R_targets)
    R0 = R_targets(r_idx);
    
    % -------- 构造回波信号 --------
    tau = 2*R0/c;
    t_eff = t - tau;
    valid = (t_eff >= 0) & (t_eff <= T);
    
    s_rx_clean = zeros(size(t));
    s_rx_clean(valid) = exp(1j*pi*u*t_eff(valid).^2);  % 去掉多普勒
    
    Pr = (Pt*Gt*Gr*lambda^2*sigma)/((4*pi)^3*R0^4*Lsys);
    s_rx_clean = sqrt(Pr)*s_rx_clean;
    
    % -------- 构造噪声用于脉压后 SNR --------
    Pn = k*T0*B_rx*NF;
    noise = sqrt(Pn/2)*(randn(size(t))+1j*randn(size(t)));
    y_n = conv(noise, conj(flip(s_tx)), 'full');
    
    % -------- 循环频偏 --------
    snr_vs_df = zeros(size(df_list));
    for k_df = 1:length(df_list)
        df_current = df_list(k_df);
        h = conj(flip(s_tx .* exp(-1j*2*pi*df_current*t))); % 匹配滤波器频偏
        rx = s_rx_clean + noise;
        rx_filt = filter(rx_lpf, 1, rx);

        y_sig = conv(rx_filt, h, 'full');
        y_n   = conv(filter(rx_lpf,1,noise), h, 'full');

%         y_sig = conv(s_rx_clean, h, 'full');
        [~, idx] = max(abs(y_sig));
        snr_vs_df(k_df) = 10*log10(abs(y_sig(idx))^2 / var(y_n));
        
    end
   %% ===== 功能级 SNR 下降量（频段重合比例模型）=====
    B_sig = B;          % 信号带宽
    B_rx_func = B;      % 功能级假设接收机带宽 = 信号带宽

    overlap_ratio = zeros(size(df_list));
    snr_drop_func = zeros(size(df_list));

    for k_df = 1:length(df_list)
        df = df_list(k_df);

        % 发射频段（相对 fc）
        Ftl = -B_sig/2 + df;
        Ftu =  B_sig/2 + df;

        % 接收机频段
        Frl = -B_rx_func/2;
        Fru =  B_rx_func/2;

        % 频段重合长度
        B_overlap = max(0, min(Ftu, Fru) - max(Ftl, Frl));
        overlap_ratio(k_df) = B_overlap / B_rx_func;
        % 功能级 SNR 下降量
        snr_drop_func(k_df) = 10*log10(B_overlap / B_rx_func + eps);
    end

    snr0 = max(snr_vs_df);
    snr_drop_sig = snr_vs_df - snr0;
    df_norm = df_list / B;
    plot(df_norm, snr_drop_sig, '.-','Color',colors(r_idx,:), ...
     'LineWidth',1.5, ...
     'DisplayName',sprintf('信号级 R = %d m', R0));

end
   valid_idx = overlap_ratio > 0;

plot(df_norm(valid_idx), snr_drop_func(valid_idx), 'k--','LineWidth',2, ...
     'DisplayName','功能级估算（部分频段重合区）');

xline(1, 'k:', 'LineWidth', 1.8, ...
      'DisplayName', '完全频段失配 (\Delta f = B)');

xlabel('归一化频偏 \Delta f / B');
ylabel('SNR 下降量 \DeltaSNR (dB)');
title('基于频段重合比例的信号级与功能级 SNR 下降对比');
legend('Location','best');
grid on;


%% ================= 多频偏 SNR vs 距离 =================
R_list = linspace(50, 4000, 100); 
df_plot_list = [0 0.5e6 1e6 2e6 12e6];  % 多个频偏
colors = lines(length(df_plot_list));

figure; hold on;
for k_df = 1:length(df_plot_list)
    df = df_plot_list(k_df);
    
    snr_sig_dist = zeros(size(R_list));
    
    for i = 1:length(R_list)
        R0 = R_list(i);
        tau = 2*R0/c;
        t_eff = t - tau;
        valid = (t_eff >= 0) & (t_eff <= T);
        
        % 回波
        s_rx_clean = zeros(size(t));
        s_rx_clean(valid) = exp(1j*pi*u*t_eff(valid).^2);
        Pr = (Pt*Gt*Gr*lambda^2*sigma)/((4*pi)^3*R0^4*Lsys);
        s_rx_clean = sqrt(Pr)*s_rx_clean;
        
        % 噪声
        Pn = k*T0*B_rx*NF;
        noise = sqrt(Pn/2)*(randn(size(t))+1j*randn(size(t)));
        
        % 匹配滤波
        h_df = conj(flip(s_tx .* exp(-1j*2*pi*df*t)));
        y_sig = conv(s_rx_clean + noise, h_df, 'full');
        y_n   = conv(noise, h_df, 'full');
        [~, idx] = max(abs(y_sig));
        snr_sig_dist(i) = 10*log10(abs(y_sig(idx))^2 / var(y_n));
    end
    
    plot(R_list, snr_sig_dist, 'Color', colors(k_df,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\Delta f = %.1f MHz', df/1e6));
end

xlabel('目标距离 R (m)');
ylabel('信号级 SNR (dB)');
title('不同频偏下信号级 SNR 随距离变化');
legend('Location','best');
grid on;


%% ================= 多频偏距离误差 vs 距离 =================
figure; hold on;
for k_df = 1:length(df_plot_list)
    df = df_plot_list(k_df);
    R_err_list = zeros(size(R_list));
    
    for i = 1:length(R_list)
        R0 = R_list(i);
        tau = 2*R0/c;
        t_eff = t - tau;
        valid = (t_eff >= 0) & (t_eff <= T);

        s_rx_clean = zeros(size(t));
        s_rx_clean(valid) = exp(1j*pi*u*t_eff(valid).^2);
        Pr = (Pt*Gt*Gr*lambda^2*sigma)/((4*pi)^3*R0^4*Lsys);
        s_rx_clean = sqrt(Pr)*s_rx_clean;

        Pn = k*T0*B_rx*NF;
        noise = sqrt(Pn/2)*(randn(size(t))+1j*randn(size(t)));

        h_df = conj(flip(s_tx .* exp(-1j*2*pi*df*t)));
        y_sig = conv(s_rx_clean + noise, h_df, 'full');
        [~, idx] = max(abs(y_sig));
        tau_hat = (idx - N)/fs;
        R_hat = c * tau_hat / 2;
        R_err_list(i) = R_hat - R0;
    end
    
    plot(R_list, abs(R_err_list), 'Color', colors(k_df,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\Delta f = %.1f MHz', df/1e6));
end

% ±2%阈值
yline(0.02*R_list(end), 'k--', 'LineWidth', 1.2, 'DisplayName', '2% 阈值');
xlabel('目标距离 R (m)');
ylabel('距离估计绝对误差 |R_{err}| (m)');
title('不同频偏下距离估计误差随距离变化');
legend('Location','best');
grid on;
