close all; clear; clc;

%% ================= LFM 信号参数 =================
B = 2e6;
T = 100e-6;
u = B / T;
fs = 4 * B;
N = round(T * fs);
t = (0:N-1).' / fs;

%% ================= 雷达与噪声参数 =================
Pt = 1e3;
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

%% ================= 距离扫描 =================
R_list = linspace(50, 6000, 100);

snr_sig  = zeros(size(R_list));
snr_func = zeros(size(R_list));
snr_coh  = zeros(size(R_list));

%% ================= 相干积分参数 =================
N_pulse = 10;  % 积分脉冲数
sigma_phi = deg2rad(10);  % 相位抖动标准差（例如 10°）

%% ================= 发射信号 & 匹配滤波 =================
s_tx = exp(1j*pi*u*t.^2) .* exp(1j*2*pi*fc_radar*t);
h = conj(flip(s_tx));  % 匹配滤波器

for i = 1:length(R_list)
    snr_sig_pre_acc = 0;
    snr_sig_post_acc = 0;
    m_pc_eq_acc = 0;
    R_err_acc = 0;

    R0 = R_list(i);
    tau = 2*R0/c;
    if tau >= T
        snr_sig(i) = NaN; snr_func(i) = NaN; snr_coh(i) = NaN;
        continue;
    end

    %% -------- 回波信号 --------
    t_eff = t - tau;
    valid = (t_eff >= 0) & (t_eff <= T);
    s_rx_clean = zeros(size(t));
    s_rx_clean(valid) = exp(1j*pi*u*t_eff(valid).^2);  

    %% -------- 雷达方程 --------
    Pr = (Pt*Gt*Gr*lambda^2*sigma)/((4*pi)^3*R0^4*Lsys);
    s_rx_clean = sqrt(Pr)*s_rx_clean;
        


    %% -------- 单脉冲信号加噪 --------
    Pn = k*T0*B_rx*NF;
    noise = sqrt(Pn/2)*(randn(size(t))+1j*randn(size(t)));
    s_rx = s_rx_clean + noise;

    %% -------- 单脉冲匹配滤波 --------
    y_sig = conv(s_rx_clean, h, 'full');
    y     = conv(s_rx,       h, 'full');
    y_n   = conv(noise,      h, 'full');

    % [~, idx] = max(abs(y_sig));
    % snr_post = abs(y_sig(idx))^2 / var(y_n);
    % snr_sig(i) = 10*log10(snr_post);
    mag2 = abs(y).^2;
    [~, idx] = max(mag2);
    main_lobe_half = round(fs/B * 2);

    k_left  = max(1, idx-main_lobe_half);
    k_right = min(length(mag2), idx+main_lobe_half);
    
    % 主瓣功率（取平均更稳定）
    P_signal = mean(mag2(k_left:k_right));
    
    % 背景区域
    mask = true(size(mag2));
    mask(k_left:k_right) = false;
    
    P_noise = mean(mag2(mask));
    
    snr_sig(i) = 10*log10(P_signal / P_noise);


    %% -------- 相干积分 --------
    y_coh = zeros(size(y));

    for p = 1:N_pulse
        % 相位抖动（模拟相干不完美）
        phi_p = sigma_phi * randn;

        noise_p = sqrt(Pn/2)*(randn(size(t))+1j*randn(size(t)));

        % 带相位误差的回波
        s_rx_p = s_rx_clean * exp(1j*phi_p) + noise_p;

        y_p = conv(s_rx_p, h, 'full');
        y_coh = y_coh + y_p;
    end
    % 
    % [~, idx_coh] = max(abs(y_coh));
    % snr_coh(i) = 10*log10(abs(y_coh(idx_coh))^2 / (var(y_n)*N_pulse));
    mag2_coh = abs(y_coh).^2;

    [~, idx_coh] = max(mag2_coh);
    
    k_left  = max(1, idx_coh-main_lobe_half);
    k_right = min(length(mag2_coh), idx_coh+main_lobe_half);
    
    P_signal_coh = mean(mag2_coh(k_left:k_right));
    
    mask = true(size(mag2_coh));
    mask(k_left:k_right) = false;
    
    P_noise_coh = mean(mag2_coh(mask));
    
    snr_coh(i) = 10*log10(P_signal_coh / P_noise_coh);
end

%% ================= 绘制 SNR 随距离变化 =================
figure;
plot(R_list, snr_sig,'r.-','LineWidth',1.5); hold on;
plot(R_list, snr_coh,'b.-','LineWidth',1.5);
xlabel('距离 (m)');
ylabel('SNR (dB)');
legend('单脉冲信号级','相干积分 SNR','Location','best');
title(sprintf('SNR 随距离变化 (N_{pulse}=%d)', N_pulse));
grid on;

%% ================= 绘制相干积分增益对比理论 =================
figure;
snr_gain = snr_coh - snr_sig;   % 相干积分增益
plot(R_list, snr_gain,'b.-','LineWidth',1.5); hold on;
plot(R_list, 10*log10(N_pulse)*ones(size(R_list)),'r--','LineWidth',1.5);
xlabel('距离 (m)');
ylabel('相干积分 SNR 提升 (dB)');
legend('仿真相干积分增益','理论 10log_{10}(N_{pulse})','Location','best');
title(sprintf('相干积分 SNR 提升随距离变化 (N_{pulse}=%d)', N_pulse));
grid on;

%% ================= 固定距离下，相干增益 vs Np =================
R_test = 500;                 % 固定距离 500 m
Np_list = 1:30;               % 扫描脉冲数
snr_gain_sim = zeros(size(Np_list));

%% -------- 重新生成该距离下的回波 --------
tau = 2*R_test/c;
t_eff = t - tau;
valid = (t_eff >= 0) & (t_eff <= T);

s_rx_clean = zeros(size(t));
s_rx_clean(valid) = exp(1j*pi*u*t_eff(valid).^2);

Pr = (Pt*Gt*Gr*lambda^2*sigma)/((4*pi)^3*R_test^4*Lsys);
s_rx_clean = sqrt(Pr)*s_rx_clean;

Pn = k*T0*B_rx*NF;

%% -------- 单脉冲基准 SNR（信号级）--------
noise0 = sqrt(Pn/2)*(randn(size(t))+1j*randn(size(t)));
y_sig0 = conv(s_rx_clean, h, 'full');
y_n0   = conv(noise0,     h, 'full');

% [~, idx0] = max(abs(y_sig0));
% snr_single = abs(y_sig0(idx0))^2 / var(y_n0);
%% -------- 单脉冲基准 SNR（主瓣/背景定义）--------

% 单脉冲输出（含噪）
noise0 = sqrt(Pn/2)*(randn(size(t))+1j*randn(size(t)));
s_rx0  = s_rx_clean + noise0;

y0 = conv(s_rx0, h, 'full');

% 功率序列
mag2 = abs(y0).^2;

% 找峰值
[~, idx0] = max(mag2);

% ===== 主瓣宽度 =====
main_lobe_half = round(fs/B * 2);   % 2倍 1/B 更稳健

k_left  = max(1, idx0-main_lobe_half);
k_right = min(length(mag2), idx0+main_lobe_half);

% ===== 主瓣功率 =====
P_signal = mean(mag2(k_left:k_right));

% ===== 背景功率 =====
mask = true(size(mag2));
mask(k_left:k_right) = false;

P_noise = mean(mag2(mask));

% ===== 单脉冲 SNR =====
snr_single = P_signal / P_noise;


%% -------- 扫描 Np --------
for ii = 1:length(Np_list)
    Np = Np_list(ii);

    y_coh = zeros(size(y_sig0));

    for p = 1:Np
        phi_p = sigma_phi * randn;  % 相位抖动
        noise_p = sqrt(Pn/2)*(randn(size(t))+1j*randn(size(t)));

        s_rx_p = s_rx_clean * exp(1j*phi_p) + noise_p;
        y_p = conv(s_rx_p, h, 'full');

        y_coh = y_coh + y_p;
    end

    % [~, idx] = max(abs(y_coh));
    % snr_coh_np = abs(y_coh(idx))^2 / (var(y_n0)*Np);
    mag2_coh = abs(y_coh).^2;

    [~, idx] = max(mag2_coh);
    
    k_left  = max(1, idx-main_lobe_half);
    k_right = min(length(mag2_coh), idx+main_lobe_half);
    
    P_signal = mean(mag2_coh(k_left:k_right));
    
    mask = true(size(mag2_coh));
    mask(k_left:k_right) = false;
    
    P_noise = mean(mag2_coh(mask));
    
    snr_coh_np = P_signal / P_noise;

    snr_gain_sim(ii) = 10*log10(snr_coh_np / snr_single);
end

%% ================= 绘图 =================
figure;
plot(Np_list, snr_gain_sim,'b.-','LineWidth',1.5); hold on;
plot(Np_list, 10*log10(Np_list),'r--','LineWidth',1.5);
xlabel('相干积分脉冲数 N_p');
ylabel('相干积分增益 (dB)');
legend('信号级仿真','理论 10log_{10}(N_p)','Location','northwest');
title('R = 500 m 处相干积分增益随 N_p 变化');
grid on;

