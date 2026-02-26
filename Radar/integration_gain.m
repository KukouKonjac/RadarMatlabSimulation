close all; clear; clc;

%% ================= 基本参数 =================
B = 2e6; T = 100e-6;
u = B/T;
fs = 4*B;
N = round(T*fs);
t = (0:N-1).'/fs;

Pt = 1e3; Gt = 30; Gr = 30;
sigma = 1; Lsys = 1;

k = physconst('Boltzmann');
T0 = 290; NF = 3;
B_rx = fs;

c = physconst('LightSpeed');
fc = 10e9;
lambda = c/fc;

Pn = k*T0*B_rx*NF;

%% ================= 信号 & 匹配滤波器 =================
s_tx = exp(1j*pi*u*t.^2).*exp(1j*2*pi*fc*t);
h = conj(flip(s_tx));

%% ================= 距离扫描 =================
R_list = linspace(50,6000,100);
N_pulse = 10;
sigma_phi = 0;
N_MC = 50;

snr_sig = zeros(size(R_list));
snr_coh = zeros(size(R_list));

for i = 1:length(R_list)

    R0 = R_list(i);
    tau = 2*R0/c;
    if tau >= T
        snr_sig(i)=NaN; 
        snr_coh(i)=NaN; 
        continue;
    end

    % ===== 回波 =====
    t_eff = t - tau;
    s_clean = zeros(size(t));
    idx_valid = (t_eff>=0)&(t_eff<=T);
    s_clean(idx_valid) = exp(1j*pi*u*t_eff(idx_valid).^2);

    Pr = (Pt*Gt*Gr*lambda^2*sigma)/((4*pi)^3*R0^4*Lsys);
    s_clean = sqrt(Pr)*s_clean;

    snr_s = zeros(1,N_MC);
    snr_c = zeros(1,N_MC);

    for mc = 1:N_MC

        % ===== 单脉冲 =====
        noise = sqrt(Pn/2)*(randn(size(t))+1j*randn(size(t)));
        y = conv(s_clean+noise,h,'full');
        snr_s(mc) = calc_snr_peak(y,fs,B);

        % ===== 相干积分 =====
        y_coh = 0;
        for p = 1:N_pulse
            phi = sigma_phi*randn;
            noise_p = sqrt(Pn/2)*(randn(size(t))+1j*randn(size(t)));
            y_coh = y_coh + ...
                conv(s_clean*exp(1j*phi)+noise_p,h,'full');
        end

        snr_c(mc) = calc_snr_peak(y_coh,fs,B);

    end

    % Monte Carlo 线性域平均
    snr_sig(i) = 10*log10(mean(snr_s));
    snr_coh(i) = 10*log10(mean(snr_c));

end

%% ================= SNR vs 距离 =================
figure;
plot(R_list,snr_sig,'r.-','LineWidth',1.5); hold on;
plot(R_list,snr_coh,'b.-','LineWidth',1.5);
xlabel('距离 (m)');
ylabel('SNR (dB)');
legend('单脉冲','相干积分','Location','best');
title(sprintf('SNR vs 距离 (N_p=%d)',N_pulse));
grid on;

%% ================= 相干增益 vs 距离 =================
figure;
snr_gain = snr_coh - snr_sig;
plot(R_list,snr_gain,'b.-','LineWidth',1.5); hold on;
plot(R_list,10*log10(N_pulse)*ones(size(R_list)),'r--','LineWidth',1.5);
xlabel('距离 (m)');
ylabel('相干增益 (dB)');
legend('仿真','理论 10log_{10}(N_p)','Location','best');
title('相干增益 vs 距离');
grid on;

%% ================= 固定距离：增益 vs Np =================
R_test = 1500;
Np_list = 1:30;
N_MC = 100;

tau = 2*R_test/c;
t_eff = t - tau;
s_clean = zeros(size(t));
idx_valid = (t_eff>=0)&(t_eff<=T);
s_clean(idx_valid) = exp(1j*pi*u*t_eff(idx_valid).^2);

Pr = (Pt*Gt*Gr*lambda^2*sigma)/((4*pi)^3*R_test^4*Lsys);
s_clean = sqrt(Pr)*s_clean;

% ===== 单脉冲基准 SNR（MC平均）=====
snr_single_mc = zeros(1,N_MC);
for mc = 1:N_MC
    noise = sqrt(Pn/2)*(randn(size(t))+1j*randn(size(t)));
    y0 = conv(s_clean+noise,h,'full');
    snr_single_mc(mc) = calc_snr_peak(y0,fs,B);
end
snr_single = mean(snr_single_mc);

snr_gain_sim = zeros(size(Np_list));

for ii = 1:length(Np_list)

    Np = Np_list(ii);
    snr_tmp = zeros(1,N_MC);

    for mc = 1:N_MC
        y_coh = 0;
        for p = 1:Np
            phi = sigma_phi*randn;
            noise_p = sqrt(Pn/2)*(randn(size(t))+1j*randn(size(t)));
            y_coh = y_coh + ...
                conv(s_clean*exp(1j*phi)+noise_p,h,'full');
        end

        snr_tmp(mc) = calc_snr_peak(y_coh,fs,B);
    end

    snr_gain_sim(ii) = 10*log10(mean(snr_tmp)/snr_single);

end

figure;
plot(Np_list,snr_gain_sim,'b.-','LineWidth',1.5); hold on;
plot(Np_list,10*log10(Np_list),'r--','LineWidth',1.5);
xlabel('N_p');
ylabel('相干增益 (dB)');
legend('仿真','理论 10log_{10}(N_p)','Location','northwest');
title('R = 1500 m 相干增益 vs N_p');
grid on;

%% ================= 峰值/背景方差 SNR 计算函数 =================
function snr_val = calc_snr_peak(y,fs,B)

    y_power = abs(y).^2;
    % 主峰
    [P_peak, k_obs] = max(y_power);

    % 主瓣保护窗口（≈5个压缩宽度）
    win_size = round(fs/B*5);

    mask = true(size(y_power));
    mask(max(1,k_obs-win_size):min(length(y_power),k_obs+win_size)) = false;

    % 背景噪声方差（复数域）
    P_background = var(y(mask));

    % SNR（线性）
    snr_val = P_peak / P_background;

end