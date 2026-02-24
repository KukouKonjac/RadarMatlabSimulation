close all; clear; clc;

%% ================= LFM 信号参数 =================
B = 20e6;        % 带宽
T = 100e-6;      % 脉冲宽度
u = B / T;       % 调频率
fs = 4 * B;      % 采样率
N = round(T * fs);
t = (0:N-1).' / fs;

%% ================= 雷达参数 =================
Pt = 1e3; Gt = 30; Gr = 30; sigma = 1; Lsys = 1;
k = physconst('Boltzmann'); T0 = 290; NF = 3; B_rx = fs;

c = physconst('LightSpeed'); 
fc_radar = 10e9; 
lambda = c / fc_radar;

%% ================= 距离扫描 =================
R_list = linspace(50, 1e4, 300);

snr_sig        = zeros(size(R_list));
snr_sig_rpj    = zeros(size(R_list));
snr_sig_cpj    = zeros(size(R_list));
R_err          = zeros(size(R_list));
R_err_rpj      = zeros(size(R_list));
R_err_cpj      = zeros(size(R_list));

%% ================= 发射信号 & 匹配滤波 =================
s_tx = exp(1j*pi*u*t.^2);
h = conj(flip(s_tx));  % 匹配滤波器

%% ================= 干扰机参数 =================
P_j = 2e-10;             % 干扰总功率
delta_j = 0;            % 干扰中心频率
B_j = 20e6;             % 干扰带宽
tau_j = T/5;            % RPJ 脉冲宽度
num_pulses = 5;         % RPJ 脉冲数量

% %% ---------------- 随机脉冲干扰（RPJ） ----------------
% s_rpj = zeros(size(t));
% for n = 1:num_pulses
%     start_idx = randi([1 N-round(tau_j*fs)]);
%     pulse_idx = start_idx : start_idx + round(tau_j*fs) - 1;
%     % 分配总功率到每个脉冲
% %     P_on = P_j*T/(num_pulses*tau_j);
%     P_on = P_j/(num_pulses);
%     s_rpj(pulse_idx) = sqrt(P_on) * exp(1j*2*pi*(delta_j + (rand-0.5)*B_j)*t(pulse_idx));
% end
% 
% %% ---------------- 覆盖脉冲干扰（CPJ） ----------------
% s_cpj = sqrt(P_j) * exp(1j*2*pi*delta_j*t);  % 覆盖整个脉冲
% % 若需要带宽可加 LFM 调制：
% % s_cpj = sqrt(P_j) * exp(1j*pi*B_j/T*t.^2 + 1j*2*pi*fc_j*t);

MC = 10;   % 建议 100~500，RPJ 至少 200

snr_sig_mc     = zeros(MC,1);
snr_sig_rpj_mc = zeros(MC,1);
snr_sig_cpj_mc = zeros(MC,1);

R_err_mc       = zeros(MC,1);
R_err_rpj_mc   = zeros(MC,1);
R_err_cpj_mc   = zeros(MC,1);

%% ================= 主循环 =================
for i = 1:length(R_list)
    R0 = R_list(i);
    tau = 2*R0/c;
    if tau >= T
        snr_sig(i) = NaN; snr_sig_rpj(i) = NaN; snr_sig_cpj(i) = NaN;
        R_err(i) = NaN; R_err_rpj(i) = NaN; R_err_cpj(i) = NaN;
        continue;
    end

    for mc = 1:MC

        %% ===== 每次 MC 都重新生成噪声 =====
        Pn = k*T0*B_rx*NF;
        noise = sqrt(Pn/2)*(randn(size(t)) + 1j*randn(size(t)));

        %% ===== 每次 MC 都重新生成 RPJ =====
        s_rpj = zeros(size(t));
        for n = 1:num_pulses
            start_idx = randi([1 N-round(tau_j*fs)]);
            pulse_idx = start_idx : start_idx + round(tau_j*fs) - 1;
%             P_on = P_j/num_pulses;
            P_on = P_j*T/(num_pulses*tau_j);
            s_rpj(pulse_idx) = sqrt(P_on) * ...
                exp(1j*2*pi*(delta_j + (rand-0.5)*B_j)*t(pulse_idx));
        end

        %% ===== CPJ（可选 MC，但保持一致） =====
        s_cpj = sqrt(P_j) * exp(1j*2*pi*delta_j*t);

        %% ===== 信号回波（可放在 MC 外，但这样写最清晰） =====
        t_eff = t - tau;
        valid = (t_eff >= 0) & (t_eff <= T);
        s_rx_clean = zeros(size(t));
        s_rx_clean(valid) = exp(1j*pi*u*t_eff(valid).^2);
        Pr = (Pt*Gt*Gr*lambda^2*sigma)/((4*pi)^3*R0^4*Lsys);
        s_rx_clean = sqrt(Pr)*s_rx_clean;

        % %% ===== 干净 SNR =====
        % y_sig_clean = conv(s_rx_clean, h, 'full');
        % y_n = conv(noise, h, 'full');
        % [~, idx] = max(abs(y_sig_clean));
        % snr_sig_mc(mc) = 10*log10(abs(y_sig_clean(idx))^2 / var(y_n));
        % 
        % y_total_clean = conv(s_rx_clean + noise, h, 'full');
        % [val, k0] = max(abs(y_total_clean));
        % 
        % %% ===== 干净测距 =====
        % y_sig = conv(s_rx_clean + noise, h, 'full');
        % [~, k0] = max(abs(y_sig));

        %% ===== 1. 情况1：干净信号 + 热噪声 =====
        s_rx_total_clean = s_rx_clean + noise;
        y_total_clean = conv(s_rx_total_clean, h, 'full');
        [val, k0] = max(abs(y_total_clean)); % 寻找探测到的最高峰

        % 真实雷达视角：排除信号主瓣区域，剩下的全是"广义噪声"
        mask = true(size(y_total_clean));
        % 窗口大小可以设为脉压压缩后的理论宽度（fs/B 的几倍）
        win_size = round(fs/B * 5); 
        mask(max(1, k0-win_size):min(length(y_total_clean), k0+win_size)) = false;

        noise_floor = var(y_total_clean(mask));
        snr_sig_mc(mc) = 10*log10(val^2 / noise_floor);
        
        % 测距：基于观测到的最高峰
        tau_est = (k0 - (N-1))/fs;
        R_err_mc(mc) = c*tau_est/2 - R0;

%         %% ===== RPJ =====
%         y_rpj = conv(s_rx_clean + s_rpj, h, 'full');
%         y_n_rpj = conv(noise, h, 'full');
% %         snr_sig_rpj_mc(mc) = 10*log10(abs(y_sig_clean(idx))^2 / var(y_n_rpj));
% 
%         [~, k0_rpj] = max(abs(y_rpj));
%         tau_est_rpj = (k0_rpj - (N-1))/fs;
%         R_err_rpj_mc(mc) = c*tau_est_rpj/2 - R0;
%         snr_sig_rpj_mc(mc) = 10*log10(abs(y_rpj(k0_rpj))^2 / var(y_n_rpj));


        %% ===== 2. 情况2：RPJ 干扰 (信号 + 噪声 + 干扰) =====
        s_rx_total_rpj = s_rx_clean + s_rpj + noise;
        y_total_rpj = conv(s_rx_total_rpj, h, 'full');
        [val_rpj, k0_rpj] = max(abs(y_total_rpj)); % 寻找此时最高峰
        
        mask_rpj = true(size(y_total_rpj));
        mask_rpj(max(1, k0_rpj-win_size):min(length(y_total_rpj), k0_rpj+win_size)) = false;
        
        noise_floor_rpj = var(y_total_rpj(mask_rpj));
        snr_sig_rpj_mc(mc) = 10*log10(val_rpj^2 / noise_floor_rpj);
        
        % 测距：如果干扰峰值高于真实目标，雷达就会锁定错误距离
        tau_est_rpj = (k0_rpj - (N-1))/fs;
        R_err_rpj_mc(mc) = c*tau_est_rpj/2 - R0;

%         %% ===== CPJ =====
%         y_cpj = conv(s_rx_clean + s_cpj, h, 'full');
%         y_n_cpj = conv(noise, h, 'full');
% %         snr_sig_cpj_mc(mc) = 10*log10(abs(y_sig_clean(idx))^2 / var(y_n_cpj));
% 
%         [~, k0_cpj] = max(abs(y_cpj));
%         tau_est_cpj = (k0_cpj - (N-1))/fs;
%         R_err_cpj_mc(mc) = c*tau_est_cpj/2 - R0;
%         snr_sig_cpj_mc(mc) = 10*log10(abs(y_cpj(k0_cpj))^2 / var(y_n_cpj));

        %% ===== 3. 情况3：CPJ 干扰 (信号 + 噪声 + 干扰) =====
        s_rx_total_cpj = s_rx_clean + s_cpj + noise;
        y_total_cpj = conv(s_rx_total_cpj, h, 'full');
        [val_cpj, k0_cpj] = max(abs(y_total_cpj));
        
        mask_cpj = true(size(y_total_cpj));
        mask_cpj(max(1, k0_cpj-win_size):min(length(y_total_cpj), k0_cpj+win_size)) = false;
        
        noise_floor_cpj = var(y_total_cpj(mask_cpj));
        snr_sig_cpj_mc(mc) = 10*log10(val_cpj^2 / noise_floor_cpj);
        
        tau_est_cpj = (k0_cpj - (N-1))/fs;
        R_err_cpj_mc(mc) = c*tau_est_cpj/2 - R0;

    end

    %% ===== MC 平均（这是功能级要用的量） =====
    snr_sig(i)     = mean(snr_sig_mc);
    snr_sig_rpj(i) = mean(snr_sig_rpj_mc);
    snr_sig_cpj(i) = mean(snr_sig_cpj_mc);

    R_err(i)       = mean(R_err_mc);
    R_err_rpj(i)   = mean(R_err_rpj_mc);
    R_err_cpj(i)   = mean(R_err_cpj_mc);
end


%% ================= 图 1：SNR 对比 =================
figure; hold on;
plot(R_list, snr_sig, 'r-', 'LineWidth', 1.5);
plot(R_list, snr_sig_rpj, 'm--', 'LineWidth', 1.5);
plot(R_list, snr_sig_cpj, 'b-.', 'LineWidth', 1.5);
xlabel('距离 (m)');
ylabel('SNR (dB)');
legend({'SNR_{post-PC}', 'SNR_{post-PC, RPJ}', 'SNR_{post-PC, CPJ}'}, 'Location','best');
title('匹配滤波后 SNR 对比（干扰前 / RPJ / CPJ）');
grid on;

%% ================= 图 2：距离误差对比=================
figure; hold on;
plot(R_list, R_err, 'g-', 'LineWidth', 1.5);
plot(R_list, R_err_rpj, 'm', 'LineWidth', 1.5);
plot(R_list, R_err_cpj, 'b-.', 'LineWidth', 1.5);
xlabel('距离 (m)');
ylabel('距离估计误差 (m)');
legend({'R_{err, clean}', 'R_{err, RPJ}', 'R_{err, CPJ}'}, 'Location','best');
title('距离估计误差随距离变化（干扰前 / RPJ / CPJ）');
grid on;

%% ================= 匹配滤波后的时域波形可视化 =================
figure;

% 计算一个中间距离点，用于显示
idx_mid = round(140);
R0 = R_list(idx_mid);
tau = 2*R0/c;
t_eff = t - tau;
valid = (t_eff >= 0) & (t_eff <= T);
s_rx_clean = zeros(size(t));
s_rx_clean(valid) = exp(1j*pi*u*t_eff(valid).^2);  
Pr = (Pt*Gt*Gr*lambda^2*sigma)/((4*pi)^3*R0^4*Lsys);
s_rx_clean = sqrt(Pr)*s_rx_clean;

Pn = k*T0*B_rx*NF;
noise = sqrt(Pn/2)*(randn(size(t)) + 1j*randn(size(t)));

% 匹配滤波
y_clean = conv(s_rx_clean + noise, h, 'full');
y_rpj   = conv(s_rx_clean + noise + s_rpj, h, 'full');
y_cpj   = conv(s_rx_clean + noise + s_cpj, h, 'full');

t_y = (0:length(y_clean)-1)/fs;  % 匹配滤波后时间轴

% ------------------- subplot 1: 信号原始 -------------------
subplot(3,1,1);
plot(t_y*1e6, abs(y_clean), 'b', 'LineWidth', 1.2);
xlabel('时间 (μs)'); ylabel('|y(t)|');
title('匹配滤波后：干净信号');
grid on;

% ------------------- subplot 2: RPJ -------------------
subplot(3,1,2);
plot(t_y*1e6, abs(y_rpj), 'm', 'LineWidth', 1.2);
xlabel('时间 (μs)'); ylabel('|y(t)|');
title('匹配滤波后：随机脉冲干扰（RPJ）');
grid on;

% ------------------- subplot 3: CPJ -------------------
subplot(3,1,3);
plot(t_y*1e6, abs(y_cpj), 'b-.', 'LineWidth', 1.2);
xlabel('时间 (μs)'); ylabel('|y(t)|');
title('匹配滤波后：覆盖脉冲干扰（CPJ）');
grid on;

sgtitle('匹配滤波后的时域波形对比');

