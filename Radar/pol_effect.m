close all; clear; clc;

rng(20260130);

%% ================= LFM 信号参数 =================

B = 20e6;
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

% 雷达发射极化 (Jones Vector): 假设为水平极化
pol_tx = [1; 0]; 
% 雷达接收天线极化: 假设与发射一致
pol_rx = [1; 0];

%% ================= POL 干扰机参数 =================

PJ = 1e-10; % 到达雷达接收端的干扰功率（W）
f_pol = 5; % 极化调制频率
pol_depth = 0.8; % 极化调制深度（0~1）
theta_J = (pol_depth * pi/2) * cos(2*pi*f_pol*t);

%% ================= 距离扫描 =================

R_list = linspace(50,6000,300);
snr_sig = zeros(size(R_list));
snr_sig_pol = zeros(size(R_list));
R_err = zeros(size(R_list));
R_err_pol = zeros(size(R_list));

%% ================= 发射信号 & 匹配滤波 =================

s_tx = exp(1j*pi*u*t.^2);
h = conj(flip(s_tx));
for i = 1:length(R_list)
    R0 = R_list(i);
    tau = 2*R0/c;
    if tau >= T
    snr_sig(i) = NaN;
    snr_sig_pol(i) = NaN;
    R_err(i) = NaN;
    R_err_pol(i) = NaN;
    continue;
    
    end
    
    %% ---------- 目标回波 ----------
    
    t_eff = t - tau;
    valid = (t_eff >= 0) & (t_eff <= T);
    s_target = zeros(size(t));
    s_target(valid) = exp(1j*pi*u*t_eff(valid).^2);
    Pr = (Pt*Gt*Gr*lambda^2*sigma)/((4*pi)^3*R0^4*Lsys);
    s_target = sqrt(Pr) * s_target;
    % 假设目标是不改变极化的点目标，接收端收到的矢量信号为：
    s_target_vec = pol_tx * s_target.'; % [2 x N] 矩阵
    
    %% ---------- 热噪声 ----------
    
    Pn = k*T0*B_rx*NF;
    noise = sqrt(Pn/2)*(randn(size(t))+1j*randn(size(t)));
    
    %% ======================================================
    
    %% 情况 1：无 POL 干扰 (修正版)
    
    % 1. 获取接收机收到的标量目标信号 (考虑极化损失)
    % s_target_vec 是 [2 x N]，pol_rx' 是 [1 x 2]
    % 这里计算的是：发射极化 -> 目标反射 -> 接收天线过滤 后的总增益
    s_target_rec = pol_rx' * s_target_vec; 
    
    % 2. 转换回列向量以匹配噪声维度
    s_target_rec = s_target_rec.'; 
    
    % 3. 合成接收信号
    s_rx = s_target_rec + noise;
    
    % 4. 脉冲压缩处理 (注意：y_sig 也要用极化后的信号计算)
    y_sig = conv(s_target_rec, h, 'full'); 
    y_n = conv(noise, h, 'full');
    
    [~,k0] = max(abs(y_sig));
    snr_post = abs(y_sig(k0))^2 / var(y_n);
    snr_sig(i) = 10*log10(snr_post);
    
    % 5. 距离估计
    y = conv(s_rx, h, 'full');
    mag = abs(y);
    [~,k0] = max(mag);
    tau_est = (k0 - (N-1)) / fs;
    R_err(i) = c*tau_est/2 - R0;
    
    %% ======================================================
    
    %% 情况 2：加入 POL 极化调制干扰
    
    % phi_pol = 2*pi*f_pol*t;
    % alpha_pol = pol_depth*cos(phi_pol) + (1-pol_depth);
    phi_J = 2*pi*rand; % 干扰机随机初相
    % s_J = sqrt(PJ) * alpha_pol .* exp(1j*(2*pi*randn(size(t)) + phi_J));
    s_J_base = sqrt(PJ) * exp(1j*(2*pi*rand(size(t)) + phi_J));
    s_J_vec = [0.'; 1.'] .* s_J_base.'; % [2 x N]
    J_matrix_t = [cos(theta_J); sin(theta_J)];

    % 2. 计算干扰信号在该极化匹配下的幅度衰减系数 (PLF)
    % 即干扰矢量在雷达接收天线 pol_rx 上的投影
    amp_coeff_J = abs(pol_rx * J_matrix_t');
    fprintf('极化匹配幅度系数 (0~1): %.4f\n', amp_coeff_J);

    %% 3. 雷达接收机端的极化匹配（关键修改）
    % 雷达天线对准 pol_rx 方向进行接收，这会自然产生极化增益或损失
    % 接收信号 = (极化矢量 · 接收天线极化)
    s_target_rec = pol_rx' * s_target_vec; % 标量 [1 x N]
    s_J_rec = pol_rx' * s_J_vec;           % 标量 [1 x N]
    
    % s_rx_pol = s_target + s_J + noise;
    % 叠加噪声 (假设噪声在两个极化通道独立)
    s_rx_pol = s_target_rec.' + s_J_rec.' + noise;
    % fprintf("s_target_rec = %.4f\n",pol_rx' * pol_tx);

    y_sig_pol = conv(s_rx_pol, h, 'full');
    y_J_pol = conv(s_J_rec, h, 'full');
    y_n_pol = conv(noise, h, 'full');
    [~,k0] = max(abs(y_sig_pol));
    snr_post_pol = abs(y_sig_pol(k0))^2 / var(y_n_pol);
    snr_sig_pol(i) = 10*log10(snr_post_pol);
    
    % 距离估计
    y_pol = conv(s_rx_pol, h, 'full');
    mag = abs(y_pol);
    [~,k0] = max(mag);
    tau_est = (k0 - (N-1)) / fs;
    R_err_pol(i) = c*tau_est/2 - R0;
end

%% ================= 图 1：SNR 对比 =================
figure;
plot(R_list,snr_sig,'b','LineWidth',1.5); hold on;
plot(R_list,snr_sig_pol,'r--','LineWidth',1.5);
xlabel('距离 (m)');
ylabel('脉压后 SNR (dB)');
legend('无干扰','POL 极化调制干扰');
title('POL 干扰对脉压后 SNR 的影响');
grid on;
%% ================= 图 2：距离误差对比 =================
figure;
plot(R_list,R_err,'b'); hold on;
plot(R_list,R_err_pol,'r--');
xlabel('距离 (m)');
ylabel('距离估计误差 (m)');
legend('无干扰','POL 极化调制干扰');
title('POL 极化调制对距离估计稳定性的破坏');
grid on;