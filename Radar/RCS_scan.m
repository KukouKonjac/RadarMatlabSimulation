close all; clear; clc;

%% =======================================================
%% ========== F-16 工程真实 RCS 参数 =====================
%% =======================================================

sigma_ref = 1;   % 参考 RCS (m^2)

% ---- 机体基础散射 ----
sigma_nose   = 0.05 * sigma_ref;
sigma_side   = 1.20 * sigma_ref;
sigma_belly  = 0.60 * sigma_ref;

p_nose  = 6;
p_side  = 4;
p_belly = 4;

% ---- 主翼 / 进气道 ----
sigma_wing = 2.5 * sigma_ref;
yaw_lobe   = 30 * pi/180;
bw_yaw     = 10 * pi/180;

% ---- 尾喷口 ----
sigma_tail = 3.0 * sigma_ref;
bw_tail   = 12 * pi/180;

%% =======================================================
%% ================= Yaw 扫描 =============================
%% =======================================================

yaw_scan = linspace(-180,180,361) * pi/180;
sigma_yaw = zeros(size(yaw_scan));

for k = 1:length(yaw_scan)

    yaw = yaw_scan(k);

    % 等效对称角
    yaw_eff = abs(yaw);
    yaw_eff = min(yaw_eff, pi - yaw_eff);

    % 方向余弦
    nx = cos(yaw_eff);
    ny = sin(yaw_eff);
    nz = 0;

    % ---- 机体 ----
    sigma_body = ...
        sigma_nose * abs(nx)^p_nose + ...
        sigma_side * abs(ny)^p_side;

    % ---- 主翼 / 进气道 ----
    sigma_wing_lobe = sigma_wing * ( ...
        exp(-((yaw_eff - yaw_lobe)/bw_yaw)^2) + ...
        exp(-((yaw_eff + yaw_lobe)/bw_yaw)^2) );

    % ---- 尾喷口 ----
    sigma_tail_lobe = sigma_tail * ...
        exp(-((yaw_eff - pi)/bw_tail)^2);

    sigma_yaw(k) = sigma_body + sigma_wing_lobe + sigma_tail_lobe;
end

figure;
plot(yaw_scan*180/pi, 10*log10(sigma_yaw), 'LineWidth',1.6);
xlabel('Yaw (deg)');
ylabel('RCS (dBsm)');
title('F-16 Yaw–RCS（左右对称，含尾喷口）');
grid on;

%% =======================================================
%% ================= Pitch 扫描 ===========================
%% =======================================================

pitch_scan = linspace(-90,90,361) * pi/180;
sigma_pitch = zeros(size(pitch_scan));

for k = 1:length(pitch_scan)

    pitch = pitch_scan(k);
    pitch_eff = abs(pitch);

    % 方向余弦
    nx = cos(pitch_eff);
    ny = 0;
    nz = sin(pitch_eff);

    % ---- 机体 ----
    sigma_body = ...
        sigma_nose  * abs(nx)^p_nose + ...
        sigma_belly * abs(nz)^p_belly;

    % ---- 机腹 / 进气道增强 ----
    belly_lobe = 0.8 * sigma_ref * ...
        exp(-((pitch_eff - 40*pi/180)/(15*pi/180))^2);

    sigma_pitch(k) = sigma_body + belly_lobe;
end

figure;
plot(pitch_scan*180/pi, 10*log10(sigma_pitch), 'LineWidth',1.6);
xlabel('Pitch (deg)');
ylabel('RCS (dBsm)');
title('F-16 Pitch–RCS（机腹增强）');
grid on;


%% ================= 雷达 / 目标位置 =================
radar_pos  = [0; 0; 0];          % 雷达坐标 (m)
target_pos = [3000; 1000; 500];  % 目标坐标 (m)

%% ================= 目标姿态 =========================
yaw   = 40*pi/180;    % 偏航角
pitch = 15*pi/180;    % 俯仰角
roll  = 30*pi/180;    % 滚转角

v = target_pos - radar_pos;
R = norm(v);
v_hat = v / R;   % 雷达 -> 目标单位向量
Rz = [ cos(yaw) -sin(yaw) 0;
       sin(yaw)  cos(yaw) 0;
       0         0        1 ];

Ry = [ cos(pitch) 0 sin(pitch);
       0          1 0;
      -sin(pitch) 0 cos(pitch) ];

Rx = [ 1 0 0;
       0 cos(roll) -sin(roll);
       0 sin(roll)  cos(roll) ];

R_body = Rz * Ry * Rx;

v_body = R_body.' * v_hat;

nx = abs(v_body(1));
ny = abs(v_body(2));
nz = abs(v_body(3));
% ---- 机体 ----
sigma_body = ...
    sigma_nose  * nx^p_nose + ...
    sigma_side  * ny^p_side + ...
    sigma_belly * nz^p_belly;

% ---- 主翼 ----
yaw_eff = atan2(ny, nx);
yaw_eff = min(yaw_eff, pi - yaw_eff);

sigma_wing_lobe = sigma_wing * ( ...
    exp(-((yaw_eff - yaw_lobe)/bw_yaw)^2) + ...
    exp(-((yaw_eff + yaw_lobe)/bw_yaw)^2) );

% ---- 尾喷口 ----
sigma_tail_lobe = sigma_tail * ...
    exp(-((yaw_eff - pi)/bw_tail)^2);

sigma_current = sigma_body + sigma_wing_lobe + sigma_tail_lobe;

fprintf('当前姿态下等效 RCS = %.2f m^2 (%.2f dBsm)\n', ...
        sigma_current, 10*log10(sigma_current));
figure; hold on; axis equal; grid on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('雷达–目标相对位置与目标姿态');

% 雷达
plot3(radar_pos(1), radar_pos(2), radar_pos(3), ...
      'ks', 'MarkerFaceColor','k', 'MarkerSize',8);

% 目标
plot3(target_pos(1), target_pos(2), target_pos(3), ...
      'ro', 'MarkerFaceColor','r', 'MarkerSize',8);

% 视线
plot3([radar_pos(1) target_pos(1)], ...
      [radar_pos(2) target_pos(2)], ...
      [radar_pos(3) target_pos(3)], ...
      'k--','LineWidth',1.2);

% 目标体坐标轴
L = 300;
axes_body = R_body * eye(3);

quiver3(target_pos(1),target_pos(2),target_pos(3), ...
        L*axes_body(1,1),L*axes_body(2,1),L*axes_body(3,1), ...
        'r','LineWidth',2);

quiver3(target_pos(1),target_pos(2),target_pos(3), ...
        L*axes_body(1,2),L*axes_body(2,2),L*axes_body(3,2), ...
        'g','LineWidth',2);

quiver3(target_pos(1),target_pos(2),target_pos(3), ...
        L*axes_body(1,3),L*axes_body(2,3),L*axes_body(3,3), ...
        'b','LineWidth',2);

legend('Radar','Target','LOS','Body X','Body Y','Body Z');
view(35,25);
