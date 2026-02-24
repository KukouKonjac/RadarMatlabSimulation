close all; clear; clc;

%% ================= LFM 信号参数 =================
B = 20e6;
T = 100e-6;
u = B / T;
fs = 4 * B;
N = round(T * fs);
t = (0:N-1).' / fs;

%% ================= 雷达与噪声参数 =================
Pt = 1e10;
Gt = 30; Gr = 30;
sigma0 = 1;
Lsys = 1;

k = physconst('Boltzmann');
T0 = 290;
NF = 3;
B_rx = fs;

c = physconst('LightSpeed');
fc = 10e9;
lambda = c / fc;

%% ================= 目标运动与姿态参数 =================
v = 5000;
fd = 2*v/lambda;

target_pos = [4000; 1500; 800];
roll  = 10 * pi/180;
pitch = 75 * pi/180;
yaw   = 50 * pi/180;

p_az = 4;
p_el = 3;

%% ================= 雷达视线角 =================
x = target_pos(1);
y = target_pos(2);
z = target_pos(3);

az_los = atan2(y, x);
el_los = atan2(z, sqrt(x^2 + y^2));

d_az = az_los - yaw;
d_el = el_los - pitch;

% ===== roll 散射调制（飞机特有）=====
A_roll = 0.6;      % roll 引起的RCS起伏幅度
m_roll = 4;        % lobes 锐度
phi0   = 0;        % 主翼平面对准雷达方向

G_roll = 1 + A_roll * abs(sin(roll - phi0))^m_roll;

sigma_aspect = sigma0 * ...
    abs(sin(d_az))^p_az * ...
    abs(sin(d_el))^p_el * ...
    G_roll;


%% ================= 距离扫描 =================
R_list = linspace(500, 4000, 100);

snr_sig_const  = zeros(size(R_list));
snr_sig_aspect = zeros(size(R_list));

%% ================= Monte Carlo 次数 =================
N_mc = 20;   % 你可以调到 500 / 1000

%% ================= 发射信号 & 匹配滤波 =================
s_tx = exp(1j*pi*u*t.^2);
h = conj(flip( s_tx .* exp(1j*2*pi*fd*t) ));

Pn = k*T0*B_rx*NF;

%% ================= 主循环（信号级 + 统计） =================
for i = 1:length(R_list)

    R = R_list(i);
    tau = 2*R/c;
    if tau >= T
        snr_sig_const(i)  = NaN;
        snr_sig_aspect(i)= NaN;
        continue;
    end

    %% -------- 回波基波形 --------
    t_eff = t - tau;
    valid = (t_eff >= 0) & (t_eff <= T);

    s_rx_base = zeros(size(t));
    s_rx_base(valid) = exp(1j*pi*u*t_eff(valid).^2) ...
                     .* exp(1j*2*pi*fd*t(valid));

    %% -------- 雷达方程 --------
    Pr_const  = (Pt*Gt*Gr*lambda^2*sigma0) ...
              /((4*pi)^3*R^4*Lsys);

    Pr_aspect = (Pt*Gt*Gr*lambda^2*sigma_aspect) ...
              /((4*pi)^3*R^4*Lsys);

    s_rx_clean_const  = sqrt(Pr_const)  * s_rx_base;
    s_rx_clean_aspect = sqrt(Pr_aspect) * s_rx_base;

    %% -------- Monte Carlo 累积 --------
    snr_c_acc = 0;
    snr_a_acc = 0;

    for k_mc = 1:N_mc

        noise = sqrt(Pn/2)*(randn(size(t))+1j*randn(size(t)));

        y_c = conv(s_rx_clean_const  + noise, h, 'full');
        y_a = conv(s_rx_clean_aspect + noise, h, 'full');
        y_n = conv(noise, h, 'full');

        snr_c_acc = snr_c_acc + max(abs(y_c)).^2 / var(y_n);
        snr_a_acc = snr_a_acc + max(abs(y_a)).^2 / var(y_n);
    end

    snr_sig_const(i)   = 10*log10(snr_c_acc / N_mc);
    snr_sig_aspect(i) = 10*log10(snr_a_acc / N_mc);
end

%% ================= 图 1：SNR 对比 =================
figure;
plot(R_list, snr_sig_const, 'k','LineWidth',1.5); hold on;
plot(R_list, snr_sig_aspect,'r','LineWidth',1.5);
xlabel('距离 (m)');
ylabel('SNR (dB)');
legend('Const RCS','Aspect + Attitude RCS');
title('信号级 SNR（Monte Carlo 平滑）');
grid on;

%% ================= 图 2：三维几何与姿态示意 =================
figure; hold on; grid on; axis equal;

plot3(0,0,0,'ks','MarkerSize',10,'LineWidth',2);
plot3(x,y,z,'ro','MarkerSize',8,'LineWidth',2);
plot3([0 x],[0 y],[0 z],'k--','LineWidth',1.5);

L = 600;
vx = L*cos(yaw)*cos(pitch);
vy = L*sin(yaw)*cos(pitch);
vz = L*sin(pitch);

quiver3(x,y,z,vx,vy,vz,'r','LineWidth',2);

xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
legend('Radar','Target','Radar LOS','Target Heading');
title('雷达–目标三维相对位置与姿态');
view(35,25);
