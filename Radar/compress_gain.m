close all; clear; clc;

%% ================= 基本参数 =================
B = 20e6;  T = 100e-6;
u = B/T;
fs = 4*B;
N  = round(T*fs);
t  = (0:N-1).'/fs;

Pt = 1e3; Gt = 30; Gr = 30;
sigma = 1; Lsys = 1;

k = physconst('Boltzmann');
T0 = 290; NF = 3;
B_rx = fs;

c = physconst('LightSpeed');
fc = 10e9;
lambda = c/fc;

Pn = k*T0*B_rx*NF;

%% ================= 发射信号 & 匹配滤波 =================
s_tx = exp(1j*pi*u*t.^2);
h = conj(flip(s_tx));

%% ================= 距离扫描 =================
R_list = linspace(50,6000,300);
MC = 50;

snr_sig_pre  = zeros(size(R_list));
snr_sig_post = zeros(size(R_list));
snr_func     = zeros(size(R_list));-0
m_pc_eq      = zeros(size(R_list));
R_err        = zeros(size(R_list));

for i = 1:length(R_list)

    R0 = R_list(i);
    tau = 2*R0/c;
    if tau >= T
        snr_sig_pre(i)=NaN; snr_sig_post(i)=NaN;
        m_pc_eq(i)=NaN; R_err(i)=NaN;
        continue;
    end

    % ===== 回波信号 =====
    t_eff = t - tau;
    s_clean = zeros(size(t));
    idx = (t_eff>=0)&(t_eff<=T);

    s_clean(idx) = exp(1j*pi*u*t_eff(idx).^2);

    Pr = (Pt*Gt*Gr*lambda^2*sigma)/((4*pi)^3*R0^4*Lsys);
    s_clean = sqrt(Pr)*s_clean;

    snr_pre_acc  = 0;
    snr_post_acc = 0;
    mpc_acc      = 0;
    err_acc      = 0;

    for mc = 1:MC

        noise = sqrt(Pn/2)*(randn(size(t))+1j*randn(size(t)));
        s_rx  = s_clean + noise;

        % -------- 脉压前 SNR --------
    
        snr_pre = max(abs(s_clean).^2) / var(noise);

        % -------- 匹配滤波 --------
        y = conv(s_rx, h, 'full');
        y_power = abs(y).^2;

        % -------- 主瓣/背景 SNR --------
        [P_peak,k_obs] = max(y_power);
        win = round(fs/B*5);

        mask = true(size(y));
        mask(max(1,k_obs-win):min(end,k_obs+win)) = false;

        P_bg = var(y(mask));
        snr_post = P_peak / P_bg;

        % -------- 距离估计 --------
        
        mag = abs(y);
        if k_obs>1 && k_obs<length(mag)
            y1=mag(k_obs-1); y2=mag(k_obs); y3=mag(k_obs+1);
            delta = 0.5*(y1-y3)/(y1-2*y2+y3);
        else
            delta=0;
        end

        tau_est = ((k_obs+delta)-(N-1))/fs;
        R_est   = c*tau_est/2;

        % -------- 累加 --------
        snr_pre_acc  = snr_pre_acc  + 10*log10(snr_pre);
        snr_post_acc = snr_post_acc + 10*log10(snr_post);
        mpc_acc      = mpc_acc      + snr_post/snr_pre;
        err_acc      = err_acc      + (R_est-R0);
    end


    snr_sig_pre(i)  = snr_pre_acc/MC;
    snr_sig_post(i) = snr_post_acc/MC;
    m_pc_eq(i)      = mpc_acc/MC;
    R_err(i)        = err_acc/MC;

    snr_func(i) = 10*log10(Pr*B*T/Pn);
end

%% ================= 探测失效点 =================
idx_fail = find(abs(R_err)>20,1,'first');

if isempty(idx_fail)
    R_fail = NaN; snr_fail = NaN;
else
    R_fail = R_list(idx_fail);
    snr_fail = snr_sig_post(idx_fail);
end


%% ================= 图1：脉压增益 =================
figure;
plot(R_list,10*log10(m_pc_eq),'b','LineWidth',1.5); hold on;
yline(10*log10(B*T),'r--','LineWidth',1.5);

if ~isnan(R_fail)
    xline(R_fail,'r--','LineWidth',1.5);
    scatter(R_fail,10*log10(m_pc_eq(idx_fail)),60,'r','filled');
end

xlabel('距离 (m)');
ylabel('脉冲压缩增益 (dB)');
legend('信号级','理论 BT','失效点');
title('等效脉冲压缩增益');
grid on;


%% ================= 图2：SNR 对比 =================
figure; hold on;
plot(R_list,snr_sig_pre,'g--','LineWidth',1.2);
plot(R_list,snr_sig_post,'r','LineWidth',1.5);
plot(R_list,snr_func,'k--','LineWidth',1.2);

if ~isnan(R_fail)
    xline(R_fail,'r--','LineWidth',1.5);
    scatter(R_fail,snr_fail,60,'r','filled');
end

xlabel('距离 (m)');
ylabel('SNR (dB)');
legend('信号级 pre','信号级 post','功能级','失效点');
title('SNR 对比');
grid on;


%% ================= 图3：距离误差 =================
figure;
plot(R_list,R_err,'b','LineWidth',1.5); hold on;
if ~isnan(R_fail)
    xline(R_fail,'r--','LineWidth',1.5);
end
xlabel('距离 (m)');
ylabel('距离误差 (m)');
title('距离估计误差');
grid on;

%% ================= 图4：匹配滤波前后波形 =================
figure;
subplot(2,1,1);
plot(t*1e6,abs(s_clean),'b');
xlabel('时间 (μs)'); ylabel('|s_{rx}|');

title('匹配滤波前'); grid on;

subplot(2,1,2);
y_sig = conv(s_clean,h,'full');
plot((0:length(y_sig)-1)/fs*1e6,abs(y_sig),'r');
xlabel('时间 (μs)'); ylabel('|y|');
title('匹配滤波后'); grid on;

%% ================= 图5：不同距离主瓣示意 =================
figure; hold on;
R_plot=[200 1000 2000 3000];
cols=['r','b','g','m'];

for j=1:length(R_plot)
    R0=R_plot(j);
    tau=2*R0/c;
    t_eff=t-tau;

    s=zeros(size(t));
    idx=(t_eff>=0)&(t_eff<=T);
    s(idx)=exp(1j*pi*u*t_eff(idx).^2);
    Pr=(Pt*Gt*Gr*lambda^2*sigma)/((4*pi)^3*R0^4*Lsys);
    s=sqrt(Pr)*s;

    y=conv(s + sqrt(Pn/2)*(randn(size(t))+1j*randn(size(t))),h,'full');
    mag=abs(y);
    t_y=(0:length(y)-1)/fs;

    [~,k_obs]=max(mag);
    win=round(fs/B*5);

    plot(t_y*1e6,mag,cols(j),'LineWidth',1.2);
    xline(t_y(max(1,k_obs-win))*1e6,'--',cols(j));
    xline(t_y(min(end,k_obs+win))*1e6,'--',cols(j));
end


xlabel('时间 (μs)');
ylabel('|y(t)|');
title('匹配滤波后主瓣范围');
grid on;






















































+-legend('200m','1000m','2000m','3000m');