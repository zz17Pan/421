function rss = calculateRSS(tx_subarrays, rx_subarrays, phi_t, theta_t, phi_r, theta_r, lambda, d, N_sub)
% calculateRSS - 计算基于HSPM的物理接收信号强度 Prx = | sum_l sum_k f_l^H H_kl w_k |^2
%
% 输入参数:
%   tx_subarrays - 发射端子阵列中心位置 [NumActiveSub x 3] (仅使用索引1..NumActiveSub 对应物理 2,3,4)
%   rx_subarrays - 接收端子阵列中心位置 [NumActiveSub x 3] (仅使用索引1..NumActiveSub 对应物理 2,3,4)
%   phi_t   - 发射端估计方位角 (phi, Azimuth) (度) - 用于生成 w_k
%   theta_t - 发射端估计极角 (theta, Polar Angle from Z+) (度) - 用于生成 w_k
%   phi_r   - 接收端估计方位角 (phi, Azimuth) (度) - 用于生成 f_l
%   theta_r - 接收端估计极角 (theta, Polar Angle from Z+) (度) - 用于生成 f_l
%   lambda  - 波长 (m)
%   d       - 子阵内天线间距 (m)
%   N_sub   - 每个子阵列中的天线数量 (必须为平方数, e.g., 16 for 4x4)
%
% 输出参数:
%   rss     - 物理接收信号功率 Prx (线性值)
%
% 坐标系约定:
%   使用标准的数学球坐标系:
%   theta: 极角 (与Z+轴夹角, 0 <= theta <= 180)
%   phi: 方位角 (XY平面内与X+轴夹角, -180 <= phi < 180)
%   方向向量 u = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
%
% 模型:
%   w_k 基于 (phi_t, theta_t) 生成，满足恒模 |w_k(m,n)| = 1/2
%   f_l 基于 (phi_r, theta_r) 生成，满足恒模 |f_l(m,n)| = 1/2
%   H_kl = alpha_kl * exp(-j*delta_Phi_kl) * a_rx_true * (a_tx_true)^H
%   alpha_kl = lambda / (4*pi*D_kl) (路径增益，假设 LOS)
%   a_tx_true, a_rx_true 基于实际几何计算的真实 AoD/AoA
%   delta_Phi_kl 基于实际子阵中心距离 D_kl

% --- 基本参数 --- 
k0 = 2*pi/lambda;
k_f = 0.05; % <--- 新增: 分子吸收系数 (Nepers/m) @ 300GHz (示例值, 应根据实际频率和环境确定)
n = sqrt(N_sub); % 子阵边长
if round(n) ~= n
    error('N_sub 必须是平方数');
end
num_active_sub = size(tx_subarrays, 1); % 波束对准使用的子阵数 (应为3)
if size(rx_subarrays, 1) ~= num_active_sub
    error('发射和接收的子阵数量必须相同');
end

% --- 硬编码反射参数 (简单地面反射模型) ---
z_ground = -1.5; % 地面高度 (米)
Gamma = 0.1 * exp(1j * pi * 0.6); % <--- 修改: 更真实的复数反射系数 (幅度和相位)
enable_nlos = true; % <--- 修改: 重新启用几何 NLOS 计算

% --- 1. 生成波束赋形/合并向量 (基于估计角度) --- 
azimuth_t_rad = deg2rad(phi_t);    % 发射端方位角
polar_t_rad = deg2rad(theta_t);    % 发射端极角
azimuth_r_rad = deg2rad(phi_r);    % 接收端方位角
polar_r_rad = deg2rad(theta_r);    % 接收端极角

% 计算估计方向的阵列响应向量 (未归一化，包含相位信息)
a_tx_resp_est = compute_array_response(azimuth_t_rad, polar_t_rad, N_sub, d, k0, false); % (N_sub x 1)
a_rx_resp_est = compute_array_response(azimuth_r_rad, polar_r_rad, N_sub, d, k0, false); % (N_sub x 1)

% 生成 w_k 和 f_l (假设所有子阵使用相同的估计角度)
% Tx: 满足恒模约束 |w_k(i)| = 1/2 
W = zeros(N_sub, num_active_sub); % 每列是一个子阵的 w_k
w_element_amp = 1/2;
phase_w = angle(a_tx_resp_est); % (N_sub x 1)
w_base = w_element_amp * exp(1j * phase_w);
for k_idx = 1:num_active_sub
    W(:, k_idx) = w_base;
end

% Rx: 满足总功率约束 sum ||f_l||^2 = 1
% 假设功率均分给 num_active_sub 个子阵, 则 ||f_l||^2 = 1/num_active_sub
% f_l = base_vector / sqrt(norm(base_vector)^2 * num_active_sub)
% base_vector = exp(1j * phase_f), norm(base_vector)^2 = N_sub
F = zeros(N_sub, num_active_sub); % 每列是一个子阵的 f_l
phase_f = angle(a_rx_resp_est); % (N_sub x 1)
f_base_element_phase = exp(1j * phase_f);
f_norm_factor = 1 / sqrt(N_sub * num_active_sub); % 归一化因子使得 sum ||f_l||^2 = 1
f_base_normalized = f_norm_factor * f_base_element_phase;
for l_idx = 1:num_active_sub
    F(:, l_idx) = f_base_normalized;
end

% --- 2. 计算实际信道矩阵 H_kl (基于真实几何 - LOS + NLOS) ---

% 参考子阵 (用于计算相对相位差)
center_idx = ceil(num_active_sub / 2); % e.g., for 3 -> 2
center_tx = tx_subarrays(center_idx, :);
center_rx = rx_subarrays(center_idx, :); % 使用真实的接收端位置

% 计算 LOS 参考距离和方向 (用于相位参考)
v_ref_los = center_rx - center_tx;
D_ref_los = norm(v_ref_los);
if D_ref_los < eps
    warning('LOS参考子阵距离过近');
    D_ref_los = eps;
    u_true_tx_AoD_los = [1, 0, 0]; 
else
    u_true_tx_AoD_los = v_ref_los / D_ref_los; 
end

% 计算真实的 LOS AoD 和 AoA (所有子阵共享相同的真实角度)
% LOS AoD
true_polar_AoD_los_rad = acos(u_true_tx_AoD_los(3));          % theta from Z+
true_azimuth_AoD_los_rad = atan2(u_true_tx_AoD_los(2), u_true_tx_AoD_los(1)); % phi from X+
% LOS AoA
true_polar_AoA_los_rad = pi - true_polar_AoD_los_rad; % theta_AoA = 180 - theta_AoD
true_azimuth_AoA_los_rad = true_azimuth_AoD_los_rad + pi; % phi_AoA = phi_AoD + 180
true_azimuth_AoA_los_rad = atan2(sin(true_azimuth_AoA_los_rad), cos(true_azimuth_AoA_los_rad));

% 计算 LOS 真实阵列响应 (归一化)
a_tx_true_los = compute_array_response(true_azimuth_AoD_los_rad, true_polar_AoD_los_rad, N_sub, d, k0, true);
a_rx_true_los = compute_array_response(true_azimuth_AoA_los_rad, true_polar_AoA_los_rad, N_sub, d, k0, true);

% --- 3. 计算总接收信号 (LOS + NLOS) --- 
% y_signal = 0; % (旧的初始化，下面会重新计算)

% --- 重构计算逻辑: 在循环内合并 LOS 和 几何 NLOS ---

y_signal_total = 0; % 初始化总信号

for l_idx = 1:num_active_sub % 接收子阵索引
    rx_pos = rx_subarrays(l_idx, :); % 使用传入的真实接收端位置
    f_l = F(:, l_idx); % 当前接收子阵的合并向量 (N_sub x 1)
    
    y_signal_l = 0; % 累加到这个接收子阵的总信号

    for k_idx = 1:num_active_sub % 发射子阵索引
        tx_pos = tx_subarrays(k_idx, :);
        w_k = W(:, k_idx); % 当前发射子阵的赋形向量 (N_sub x 1)

        % --- 3.1 LOS Path Contribution --- 
        D_kl_los = norm(rx_pos - tx_pos);
        if D_kl_los < eps; D_kl_los = eps; end
        delta_phi_kl_los = k0 * (D_kl_los - D_ref_los); % 相对于LOS参考路径
        phase_term_los = exp(-1j * delta_phi_kl_los);
        % alpha_kl_los = lambda / (4*pi*D_kl_los); % (旧: 仅自由空间损耗)
        alpha_kl_los = (lambda / (4*pi*D_kl_los)) * exp(-0.5 * k_f * D_kl_los); % <--- 修改: 加入分子吸收
        H_kl_los = alpha_kl_los * phase_term_los * (a_rx_true_los * a_tx_true_los');
        y_kl_los = (f_l' * H_kl_los * w_k);

        % --- 3.2 NLOS Path Contribution (Ground Reflection - Geometric) --- 
        y_kl_nlos = 0; % 初始化NLOS贡献
        if enable_nlos && tx_pos(3) > z_ground && rx_pos(3) > z_ground % 确保都在地面之上
            % 1. Mirror image of transmitter
            tx_pos_img = [tx_pos(1), tx_pos(2), 2*z_ground - tx_pos(3)];

            % 2. Total path length
            D_kl_nlos = norm(rx_pos - tx_pos_img);
            if D_kl_nlos < eps; D_kl_nlos = eps; end

            % 3. Reflection point (for calculating angles)
            if abs(rx_pos(3) - tx_pos_img(3)) > eps
                t_refl = (z_ground - tx_pos_img(3)) / (rx_pos(3) - tx_pos_img(3));
                 if t_refl >= 0 && t_refl <= 1 % 确保物理反射点存在于 Tx_img 和 Rx 之间
                    refl_pos = tx_pos_img + t_refl * (rx_pos - tx_pos_img);

                    % 4. Calculate NLOS Angles
                    v_tx_refl = refl_pos - tx_pos;
                    v_refl_rx = rx_pos - refl_pos;
                    if norm(v_tx_refl) > eps && norm(v_refl_rx) > eps
                        u_true_tx_AoD_nlos = v_tx_refl / norm(v_tx_refl);
                        u_true_rx_AoA_dir_nlos = v_refl_rx / norm(v_refl_rx); % 指向 Rx
                        
                        % NLOS AoD angles
                        true_polar_AoD_nlos_rad = acos(u_true_tx_AoD_nlos(3));
                        true_azimuth_AoD_nlos_rad = atan2(u_true_tx_AoD_nlos(2), u_true_tx_AoD_nlos(1));

                        % NLOS AoA angles (信号到达方向的反方向的矢量)
                        u_arrival_nlos = -u_true_rx_AoA_dir_nlos; 
                        true_polar_AoA_nlos_rad = acos(u_arrival_nlos(3));
                        true_azimuth_AoA_nlos_rad = atan2(u_arrival_nlos(2), u_arrival_nlos(1));

                        % 5. Calculate NLOS Array Responses (归一化 for H matrix construction)
                        a_tx_true_nlos = compute_array_response(true_azimuth_AoD_nlos_rad, true_polar_AoD_nlos_rad, N_sub, d, k0, true); % Apply normalization
                        a_rx_true_nlos = compute_array_response(true_azimuth_AoA_nlos_rad, true_polar_AoA_nlos_rad, N_sub, d, k0, true); % Apply normalization

                        % 6. Calculate NLOS Path Gain and Phase
                        % alpha_kl_nlos = Gamma * lambda / ((4*pi*D_kl_nlos)); % (旧: 仅自由空间损耗+反射)
                        alpha_kl_nlos = Gamma * (lambda / (4*pi*D_kl_nlos)) * exp(-0.5 * k_f * D_kl_nlos); % <--- 修改: 加入分子吸收
                        delta_phi_kl_nlos = k0 * (D_kl_nlos - D_ref_los); % 相对于 *LOS* 参考路径
                        phase_term_nlos = exp(-1j * delta_phi_kl_nlos);

                        % 7. Construct NLOS Channel Matrix
                        H_kl_nlos = alpha_kl_nlos * phase_term_nlos * (a_rx_true_nlos * a_tx_true_nlos');

                        % 8. Calculate NLOS signal contribution
                        y_kl_nlos = (f_l' * H_kl_nlos * w_k);
                    end
                 end
            end
        end % End if enable_nlos

        % --- 3.3 累加 LOS 和 NLOS 贡献到子阵 l --- 
        y_signal_l = y_signal_l + y_kl_los + y_kl_nlos;

    end % End loop over k (tx_subarrays)
    
    % --- 3.4 累加来自接收子阵 l 的信号到总信号 --- 
    y_signal_total = y_signal_total + y_signal_l;
    
end % End loop over l (rx_subarrays)

% --- 4. 计算接收功率 (基于合并后的信号) ---
Prx = abs(y_signal_total)^2;

% --- 引入增益缩放因子 Beta ---
Beta = 1e12; % 新增：定义缩放因子 (可调整)
Prx_scaled = Prx * Beta; % 新增：应用缩放

% --- 返回 dB 值 (基于缩放后的功率) ---
if Prx_scaled < eps
    rss = -inf; % Or a very small dB value like -300
else
    rss = 10 * log10(Prx_scaled); % 修改：使用 Prx_scaled
end
% rss = 10 * log10(max(Prx_scaled, eps)); % Alternative to avoid -inf

end

% 辅助函数：计算阵列响应向量
function a = compute_array_response(azimuth_rad, polar_rad, N_sub, d, k0, apply_normalization)
    n = sqrt(N_sub);
    a = zeros(N_sub, 1);
    idx = 1;
    % 方向向量 u = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
    u_vec = [sin(polar_rad)*cos(azimuth_rad), sin(polar_rad)*sin(azimuth_rad), cos(polar_rad)];
    
    if apply_normalization
        norm_factor = 1/sqrt(N_sub); % Normalization factor for power calculation
    else
        norm_factor = 1; % No normalization for phase extraction
    end

    for m = 1:n % x-direction index in sub-array (对应理论 m)
        for o = 1:n % z-direction index in sub-array (对应理论 n)
            % Antenna offset relative to sub-array center (on xz plane)
            % Center at (0,0) -> (m - (n+1)/2)*d for x, (o - (n+1)/2)*d for z
            delta_r = [(m-(n+1)/2)*d, 0, (o-(n+1)/2)*d]; 
            % Phase = exp(-j * k0 * dot(delta_r, u_vec))
            phase = exp(-1j * k0 * (delta_r(1)*u_vec(1) + delta_r(3)*u_vec(3)));
            a(idx) = phase * norm_factor; % Apply normalization factor if requested
            idx = idx + 1;
        end
    end
end 