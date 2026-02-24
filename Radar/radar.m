close all; clear;

%% ================= LFM 信号参数 =================
B = 20e6;           % 带宽 20 MHz
T = 100e-6;         % 脉冲宽度 100 μs
u = B / T;          % 调频斜率
fs = 4 * B;         % 采样率
N = round(T * fs);
t = (0:N-1).' / fs;

%% ================= 雷达与噪声参数 =================
Pt = 1e3;           % 发射功率
Gt = 30;
Gr = 30;
sigma = 1;
Lsys = 1;

k = physconst('Boltzmann');
T0 = 290;
NF = 3;
B_rx = fs;

c = physconst('LightSpeed');
fc = 10e9;
lambda = c / fc;

%% ================= 目标参数 =================
v = 7000;                 % 径向速度
beta = 2 * v / c;        % 多普勒时间缩放因子
fd = 2*v/lambda;

%% ================= 距离扫描 =================
R_min = 100;
R_max = 3000;
num_ranges = 50;
R_list = linspace(R_min, R_max, num_ranges);

est_error = zeros(size(R_list));
snr_func_list = zeros(size(R_list));
snr_out_list  = zeros(size(R_list));

rng(0);

%% ================= 发射信号 =================
s_tx = exp(1j * pi * u * t.^2);
h = conj(flip(s_tx));     % 匹配滤波器
% 修改匹配滤波器为包含多普勒
h_doppler = conj(flip(s_tx)) .* exp(-1j * 2 * pi * fd * (N-1:-1:0)/fs);

%% ================= 主循环 =================
for i = 1:length(R_list)
    R0 = R_list(i);
    tau = 2 * R0 / c;

    if tau >= T
        est_error(i) = NaN;
        snr_func_list(i) = NaN;
        snr_out_list(i)  = NaN;
        continue;
    end

    %% -------- 信号级回波（连续延迟 + 多普勒时间缩放） --------
    t_rx = t;
    t_tx_eff = (1 - beta) * (t_rx - tau);

    valid_idx = (t_tx_eff >= 0) & (t_tx_eff <= T);

    s_rx_clean = zeros(size(t));
    s_rx_clean(valid_idx) = exp(1j * pi * u * t_tx_eff(valid_idx).^2);

    %% -------- 雷达方程幅度 --------
    Pr = (Pt * Gt * Gr * lambda^2 * sigma) / ...
         ((4*pi)^3 * R0^4 * Lsys);
    A = sqrt(Pr);
    s_rx_clean = A * s_rx_clean;

    %% ================= 功能级 SNR =================
    Pn = k * T0 * B_rx * NF;
    % ======脉冲压缩增益=======
    m_compress = B*T;
    % ======多普勒失配损失=====
    lambda = c / fc;
    fd = 2 * v / lambda;
    eta_d = sinc(fd * T)^2;   % MATLAB sinc(x)=sin(pi x)/(pi x)
    snr_func_list(i) = 10 * log10(Pr * m_compress / Pn) ;

    %% ================= 加热噪声 =================
    noise = sqrt(Pn/2) * ...
        (randn(size(s_rx_clean)) + 1j*randn(size(s_rx_clean)));
    s_rx = s_rx_clean + noise;

    %% ================= 匹配滤波 =================
    s_pc_sig = conv(s_rx_clean, h, 'full');  % 纯信号
    s_pc     = conv(s_rx,       h, 'full');  % 信号+噪声
    noise_pc = conv(noise,       h, 'full');  % 纯噪声

    %% ================= SNR_out 计算 =================
    [~, peak_idx] = max(abs(s_pc_sig));
    signal_peak_power = abs(s_pc_sig(peak_idx))^2;

    guard = 5;
    noise_idx = [1:peak_idx-guard, peak_idx+guard:length(noise_pc)];
    noise_var_out = var(noise_pc(noise_idx));

    snr_out_list(i) = 10 * log10(signal_peak_power / noise_var_out);

    %% ================= 距离估计 =================
    [~, peak_idx_full] = max(abs(s_pc));
    tau_est = (peak_idx_full - N) / fs;
    R_est = c * tau_est / 2;
    est_error(i) = R_est - R0;
end

%% ================= 新增：固定距离 3000m，SNR 随速度变化 =================
R_fixed = 3000;               % 固定距离
v_list = linspace(0, 10000, 100);  % 速度从 0 到 10 km/s

snr_func_v = zeros(size(v_list));
snr_out_v  = zeros(size(v_list));

% 重用之前的参数（Pt, Gt, Gr, lambda, sigma, Lsys, k, T0, NF, B_rx, fs, N, t, h）
for j = 1:length(v_list)
    v_j = v_list(j);
    
    % 计算回波延迟
    tau = 2 * R_fixed / c;
    if tau >= T
        snr_func_v(j) = NaN;
        snr_out_v(j)  = NaN;
        continue;
    end
    
    % 多普勒频移
    fd = 2 * v_j / lambda;
    
    % 发射信号时间映射（纯频移模型）
    t_tx_eff = t - tau;
    valid_idx = (t_tx_eff >= 0) & (t_tx_eff <= T);
    
    s_rx_clean = zeros(size(t));
    s_rx_clean(valid_idx) = exp(1j * pi * u * t_tx_eff(valid_idx).^2) ...
                          .* exp(1j * 2 * pi * fd * t(valid_idx));
    
    % 接收回波幅度（雷达方程）
    Pr = (Pt * Gt * Gr * lambda^2 * sigma) / ((4*pi)^3 * R_fixed^4 * Lsys);
    A = sqrt(Pr);
    s_rx_clean = A * s_rx_clean;
    
    % 功能级 SNR（含脉冲压缩增益和多普勒失配）
    Pn = k * T0 * B_rx * NF;
    m_compress = B * T;
    eta_d = sinc(fd * T)^2;   % 注意：MATLAB 的 sinc(x) = sin(pi*x)/(pi*x)
    snr_func_v(j) = 10 * log10(Pr * m_compress * eta_d / Pn);
    
    % 加噪声
    noise = sqrt(Pn/2) * (randn(size(s_rx_clean)) + 1j*randn(size(s_rx_clean)));
    s_rx = s_rx_clean + noise;
    
    % 匹配滤波
    s_pc_sig = conv(s_rx_clean, h, 'full');
    s_pc     = conv(s_rx,       h, 'full');
    noise_pc = conv(noise,      h, 'full');
    
    % 信号级 SNR（峰值功率 / 噪声方差）
    [~, peak_idx] = max(abs(s_pc_sig));
    signal_peak_power = abs(s_pc_sig(peak_idx))^2;
    
    guard = 5;
    noise_idx = [1:peak_idx-guard, peak_idx+guard:length(noise_pc)];
    noise_var_out = var(noise_pc(noise_idx));
    
    snr_out_v(j) = 10 * log10(signal_peak_power / noise_var_out);
end

%% ================= 图 1：距离误差 =================
figure;
subplot(2,1,1);
plot(R_list, est_error, 'b.-','LineWidth',1.5);
xlabel('真实距离 (m)');
ylabel('距离估计误差 (m)');
title('距离估计误差 vs 距离（基于 SNR_{out}）');
grid on;

%% ================= 图 2：SNR 对比 =================
subplot(2,1,2);
plot(R_list, snr_out_list,    'r.-','LineWidth',1.5); hold on;
plot(R_list, snr_func_list,   'k--','LineWidth',1.5);

xlabel('真实距离 (m)');
ylabel('SNR (dB)');
legend('SNR_{peak}（信号级）', ...
       'SNR_{func}（功能级）');

title('三种 SNR 定义的对比');
grid on;

%% ================= 图 3：等效增益图 =================
figure;
BT_dB = 10*log10(B*T);

plot(snr_func_list, snr_out_list, 'bo','LineWidth',1.5); hold on;
plot(snr_func_list, snr_func_list, 'r--','LineWidth',2);

xlabel('功能级 SNR_{in} (dB)');
ylabel('信号级 SNR_{out} (dB)');
legend('仿真结果','理论 SNR_{out} = SNR_{in} + 10log_{10}(BT)');
title('等效处理增益（脉冲压缩）');
grid on;

%% ================= 图 4：固定距离下 SNR 随速度变化 =================
figure;
plot(v_list, snr_out_v, 'r.-', 'LineWidth', 1.5); hold on;
plot(v_list, snr_func_v, 'k--', 'LineWidth', 1.5);
xlabel('目标径向速度 (m/s)');
ylabel('SNR (dB)');
title('固定距离 R = 3000 m 时，SNR 随速度变化（体现多普勒失配）');
legend('SNR_{out}（信号级）', 'SNR_{func}（功能级）');
grid on;
