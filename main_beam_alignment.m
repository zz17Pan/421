%% 太赫兹波束对准系统主函数
% 实现感知辅助太赫兹波束对准科研课题中的波束对准部分
% 系统模型: 固定的发射端，匀速运动的接收端，感知的目标即为接收阵列
% 阵列结构: 发射阵列与接收阵列均为4个4*4的方形阵列
%           子阵内间距为λ/2，子阵间间距为32λ
%           感知系统使用1个子阵列，波束对准使用剩余3个子阵列

clear;
close all;
clc;

%% 系统参数设置
% 基本参数
freq = 300e9;          % 工作频率：300 GHz
c = 3e8;               % 光速
lambda = c/freq;       % 波长：1 mm
d = lambda/2;          % 子阵内天线间距
d_sub = 32*lambda;     % 子阵间间距

% 阵列参数
N_sub = 16;            % 子阵内天线数 (4x4)
M_sub = 4;             % 子阵总数量
active_sub = 3;        % 波束对准使用的子阵数量（2,3,4号子阵）

% 真实位置参数 (用于生成模拟数据)
% 注意：使用标准球坐标系约定
D0 = 20;               % 发射端与接收端之间的真实距离 (m)
phi_true = 10;         % 真实方位角 phi (Azimuth) (度)
theta_true = 85;       % 真实极角 theta (Polar Angle from Z+) (度), 设为85避免90度奇点

% 生成感知系统的输出结果 (模拟值，作为波束对准系统的输入)
rng(42);  % 设置随机数种子，确保结果可重复
R0 = D0 + (rand-0.5)*1.4;           % 距离估计值，误差在±0.7m内
phi0 = phi_true + (rand-0.5)*1.6;   % 方位角估计值 (phi)，误差在±0.8°内
theta0 = theta_true + (rand-0.5)*1.6; % 极角估计值 (theta)，误差在±0.8°内
% 确保估计的极角在有效范围内 [0, 180]
theta0 = max(0, min(180, theta0));
% 标准化方位角到[-180, 180]范围
phi0 = mod(phi0 + 180, 360) - 180;

fprintf('真实值 - 距离: %.2f m, 方位角(phi): %.2f°, 极角(theta): %.2f°\n', D0, phi_true, theta_true);
fprintf('感知系统输出 - 距离: %.2f m, 方位角(phi): %.2f°, 极角(theta): %.2f°\n', R0, phi0, theta0);

%% 1. 子阵坐标建模
% 发射端子阵坐标 (xz平面, y=0)，使用2,3,4号子阵（1号子阵用于感知）
tx_subarrays = zeros(active_sub, 3);  % [x, y, z]
for k = 1:active_sub
    tx_subarrays(k, :) = [(k+1-3)*d_sub, 0, 0];  % k+1对应子阵索引2,3,4
end

%% 2. HSPM信道建模与波束成形初始化
% 初始化角度参数 (使用感知输出)
phi_t = phi0;          % 发射端方位角(phi)初始值
theta_t = theta0;      % 发射端极角(theta)初始值

% 计算接收端初始角度，使其指向发射端
% 方位角 phi_r = phi_t + 180
% 极角 theta_r = 180 - theta_t
phi_r = phi_t + 180;      % 接收端方位角(phi)初始值
phi_r = mod(phi_r + 180, 360) - 180; % 标准化到[-180, 180]范围
theta_r = 180 - theta_t;  % 接收端极角(theta)初始值
% 确保接收端极角也在[0, 180]范围内
theta_r = max(0, min(180, theta_r));

% 将初始距离设为感知估计值
R_init = R0;

%% 3. 高精度距离校准 (两阶段：窄带解模糊 + 宽带精估计)
% 使用感知系统距离值 R0 作为初始粗估计

% --- 参数定义 ---
% 宽带OFDM参数 (用于精估计)
B_fine = 4e9;             % 精估计使用带宽 4GHz
Nc = 1024;                % 子载波数量
df_fine = B_fine/Nc;      % 精估计子载波间隔
fc = freq;                % 中心频率
f_fine = fc + (-B_fine/2 : df_fine : (B_fine/2-df_fine)); % 精估计频率点

% 窄带参数 (用于粗估计解模糊)
B_coarse = 200e6;         % 粗估计使用带宽 200MHz
Nc_coarse = ceil(B_coarse / df_fine); % 保证子载波间隔一致
if mod(Nc_coarse, 2) ~= 0; Nc_coarse = Nc_coarse + 1; end % 确保偶数个点
coarse_indices = floor(Nc/2) + (-(Nc_coarse/2) : (Nc_coarse/2 - 1)); % 中心 Nc_coarse 个点
f_coarse = f_fine(coarse_indices); % 粗估计使用的频率点

% 仿真参数
snr_db = -5;              % 进一步提高信噪比，降低噪声对相位影响
snr_linear = 10^(snr_db/10);

% --- 物理测量模拟 (修正：必须基于真实距离 D0 模拟信号) ---
% 1. 模拟信号发射与接收 (在精估计的全带宽上模拟)
tx_signal = ones(size(f_fine)); % 发射信号

% 修正：模拟接收信号的相位必须基于真实距离 D0
% 双向传播相位因子：exp(-j*2π*f*2*TrueDistance/c)
rx_signal_clean = tx_signal .* exp(-1j * 2*pi*f_fine * (2*D0/c)); % 使用 D0 !!!
                                                          
% 加入噪声
noise_power = 1/snr_linear;
rx_noise = sqrt(noise_power/2) * (randn(size(f_fine)) + 1j*randn(size(f_fine)));
rx_signal = rx_signal_clean + rx_noise;

% 2. 计算测量相位 (在全带宽上)
measured_phase_fine = angle(rx_signal); % 接收信号相位 (包含噪声和模糊, 基于 D0)

% --- 阶段一：窄带粗估计 (解模糊) 使用相位一致性搜索 ---
fprintf('高精度距离估计 - 阶段一：使用相位一致性搜索解模糊 (带宽 %.1f MHz)...\n', B_coarse/1e6);

% 提取窄带部分的测量相位
measured_phase_coarse = measured_phase_fine(coarse_indices);

% 定义粗搜索范围和步长
search_range = 1.0; % 在 R0 附近 +/- 1.0 米内搜索
search_step = 0.05;  % 搜索步长 5cm
R_test_values = (R0 - search_range) : search_step : (R0 + search_range);

% 确保搜索值不为负
R_test_values = R_test_values(R_test_values > 0);

if isempty(R_test_values)
    fprintf('  警告: 搜索范围内无有效正距离，将使用 R0 作为粗估计。\n');
    R_coarse = R0;
else
    min_std_dev = inf;
    best_R_coarse = R0; % 初始化最佳值为R0
    std_devs = zeros(size(R_test_values)); % 存储每个测试距离的标准差

    % 遍历候选距离
    for i = 1:length(R_test_values)
        R_test = R_test_values(i);
        
        % 1. 使用 R_test 对窄带相位进行预补偿
        % 注意符号：测量相位 ~ -4πfR/c，所以补偿项是 +4πfR_test/c
        compensated_phase = measured_phase_coarse + (4*pi*f_coarse*R_test/c);
        
        % 2. 计算补偿后相位的循环标准差
        % 转换为复数表示，计算幅角标准差是一种常用的循环统计方法
        complex_repr = exp(1j * compensated_phase);
        mean_complex = mean(complex_repr);
        % 循环方差 V = 1 - |mean(exp(j*theta))|
        circ_variance = 1 - abs(mean_complex);
        % 循环标准差 S = sqrt(-2 * log(1-V)) = sqrt(-2 * log(abs(mean_complex)))
        % 添加小量避免 log(0)
        circ_std_dev = sqrt(-2 * log(abs(mean_complex) + 1e-12));
        std_devs(i) = circ_std_dev;
        
        % 更新最佳粗估计
        if circ_std_dev < min_std_dev
            min_std_dev = circ_std_dev;
            best_R_coarse = R_test;
        end
    end

    % 可视化搜索过程
    figure(10); clf;
    plot(R_test_values, std_devs, 'bo-');
    hold on; 
    plot(best_R_coarse, min_std_dev, 'r*', 'MarkerSize', 10);
    xlabel('测试距离 R_{test} (m)');
    ylabel('补偿后相位循环标准差 (rad)');
    title(sprintf('阶段一：相位一致性搜索 (最佳 R_{coarse} = %.4f m)', best_R_coarse));
    grid on;

    % 使用找到的最佳距离作为粗估计
    R_coarse = best_R_coarse;
    fprintf('  相位一致性搜索找到最佳粗估计 R_coarse: %.6f m (最小标准差: %.4f rad)\n', R_coarse, min_std_dev);
    
    % 合理性检查 (虽然搜索确保了正值，但还是检查一下)
    if R_coarse <= 0
       fprintf('  警告: 搜索得到的 R_coarse 非正数 (%.4f m)，回退到 R0。\n', R_coarse);
       R_coarse = R0; 
    end
end

% --- 阶段二：宽带精估计 (提精度) ---
fprintf('\n高精度距离估计 - 阶段二：宽带精估计 (使用 R_coarse=%.4f m 补偿)...\n', R_coarse);

% 1. 使用 R_coarse 对全带宽测量相位进行预补偿
% 符号修正: 测量相位是 exp(-j*2*pi*f*2*D0/c), 补偿是 exp(+j*...) 
expected_phase_coarse = -2*pi*f_fine * (2*R_coarse/c); % 基于粗估计的期望相位
residual_phase = angle(exp(1j*measured_phase_fine) .* exp(1j*-expected_phase_coarse)); % 使用复数运算避免直接加减相位

% 2. 解缠绕残差相位 (理论上应接近0，变化平缓)
unwrapped_residual = unwrap(residual_phase);

% 可视化残差相位
figure(11); clf;
subplot(2,1,1); plot(f_fine/1e9, residual_phase, 'b.'); title('宽带残差相位 (补偿后)'); ylabel('rad'); grid on;
subplot(2,1,2); plot(f_fine/1e9, unwrapped_residual, 'r.-'); title('解缠绕后的宽带残差相位'); xlabel('频率 (GHz)'); ylabel('rad'); grid on;

% 3. 线性拟合解缠绕后的残差相位，得到精细斜率
%   使用中心化和缩放改善 polyfit 的数值稳定性
[p_fine, S_fine, mu_fine] = polyfit(f_fine, unwrapped_residual, 1);
%   p_fine(1) 是基于标准化频率 xhat = (f - mu(1))/mu(2) 的斜率
%   真实斜率 = p_fine(1) / mu(2)
slope_fine = p_fine(1) / mu_fine(2); % 修正斜率计算

% 检查拟合质量
% polyval 使用 mu 来评估原始数据点的拟合值和误差
[fit_vals, delta] = polyval(p_fine, f_fine, S_fine, mu_fine); 
residual_std = std(unwrapped_residual - fit_vals); % 计算实际残差的标准差
fprintf('  宽带残差拟合斜率 (修正后): %.4e\n', slope_fine);
fprintf('  残差相位拟合标准差 (修正后): %.4f rad (越小越好)\n', residual_std);

% 4. 从精细斜率计算距离修正量 delta_R
delta_R = -slope_fine * c / (4*pi);
fprintf('  距离修正量 delta_R: %.6f m\n', delta_R);

% 5. 计算最终高精度距离估计
R_opt = R_coarse + delta_R;

% 合理性检查
if R_opt <= 0
    fprintf('  警告: 最终距离估计 %.6f m 为非正数，结果可能不可靠！回退到 R_coarse。\n', R_opt);
    R_opt = R_coarse;
end

% --- 结果输出 ---
% 计算距离估计误差 (仅用于仿真评估)
distance_error = abs(R_opt - D0);
fprintf('\n最终高精度距离估计 R_opt: %.6f m (初始感知 R0: %.4f m)\n', R_opt, R0);
fprintf('相对于真实值的误差: %.6f m (目标 < 0.05 m)\n', distance_error);

% --- 计算真实的物理接收端位置 --- 
% 注意：这里使用真实的距离 D0 和真实角度计算物理位置，因为这是静态场景
center_rx_pos_true = [D0 * sin(deg2rad(theta_true)) * cos(deg2rad(phi_true)), ...
                      D0 * sin(deg2rad(theta_true)) * sin(deg2rad(phi_true)), ...
                      D0 * cos(deg2rad(theta_true))]; 
rx_subarrays_true = zeros(active_sub, 3);  % [x, y, z]
for l = 1:active_sub
    % 阵列内部沿x轴排布，整体平移到 center_rx_pos_true
    offset = [(l+1-3)*d_sub, 0, 0]; % 假设接收子阵的内部排布与发射端相似
    rx_subarrays_true(l, :) = center_rx_pos_true + offset;
end
% ------------------------------------

%% 4. 高精度角度优化
% 设置梯度优化参数 (需要针对物理 Prx 大幅调整)
max_iter = 200;        % 增加最大迭代次数
delta_angle = 0.2;     % 梯度计算扰动步长 (度) (增大)
alpha_phi = 0.05;      % 方位角(phi)初始学习率 (减小)
alpha_theta = 0.05;    % 极角(theta)初始学习率 (减小)
converge_threshold_angle = 0.005; % 收敛角度阈值 (度)
converge_threshold_rss_abs = 1e-6; % 绝对RSS(dB)变化收敛阈值 (放宽)
min_iterations = 50;       % 最小迭代次数
grad_zero_threshold = 1e-6;  % 梯度判断阈值 (增大?)
rss_improve_threshold = 1e-5;  % RSS改善判断阈值 (dB) (放宽)
max_backtrack = 5;         % 最大回溯次数 (暂时保留)
max_stuck_count = 10;        % 最大连续无改善次数 (增加)
perturb_angle = 0.5;         % 随机扰动幅度 (度) (增大)

% 存储优化过程
phi_history = zeros(max_iter, 1);
theta_history = zeros(max_iter, 1);
rss_history = zeros(max_iter, 1);

% 当前角度估计 (方位角phi, 极角theta)
phi_current = phi0;
theta_current = theta0;

% 记录最佳值
best_phi = phi_current;
best_theta = theta_current;

% 计算当前接收端角度
phi_r = phi_current + 180;  % 方位角(phi)指向发射端
phi_r = mod(phi_r + 180, 360) - 180; % 标准化
theta_r = 180 - theta_current; % 极角(theta)反向
theta_r = max(0, min(180, theta_r)); % 限制范围

% 初始RSS计算 (使用 rx_subarrays_true)
rss_current = calculateRSS(tx_subarrays, rx_subarrays_true, phi_current, theta_current, ...
                          phi_r, theta_r, lambda, d, N_sub);

% 打印初始RSS值以便诊断 (现在是 dB)
fprintf('初始RSS值: %.3f dB\n', rss_current);

best_rss = rss_current;

% 随机重启计数
restart_count = 0;
max_restarts = 5;  % 增加重启次数

% 尝试不同初始条件
initial_points = zeros(max_restarts+1, 2);
initial_points(1,:) = [phi0, theta0];  % [方位角, 极角]

% 生成额外的初始点
for r = 1:max_restarts
    rand_azimuth = phi0 + (rand()*2-1) * 2; % ±2度内随机扰动方位角
    rand_polar = theta0 + (rand()*2-1) * 2;   % ±2度内随机扰动极角
    initial_points(r+1,:) = [
        mod(rand_azimuth + 180, 360) - 180, ... % 标准化方位角
        max(0, min(180, rand_polar)) ...     % 限制极角范围
    ];
end

% 追踪所有尝试中的最佳结果
global_best_rss = -inf;
global_best_phi = phi0;
global_best_theta = theta0;

% 主循环 - 尝试不同的初始点
for restart_idx = 1:(max_restarts+1)
    % 设置当前初始点 [方位角, 极角]
    phi_current = initial_points(restart_idx, 1);
    theta_current = initial_points(restart_idx, 2);

    % 重置学习率和状态
    current_alpha_phi = alpha_phi;
    current_alpha_theta = alpha_theta;
    consecutive_decrease = 0; % 重置连续下降计数
    grad_zero_count = 0;      % 重置梯度零计数

    % 计算初始接收端角度
    phi_r = phi_current + 180;
    phi_r = mod(phi_r + 180, 360) - 180;
    theta_r = 180 - theta_current;
    theta_r = max(0, min(180, theta_r));

    % 计算初始RSS (使用 rx_subarrays_true)
    rss_current = calculateRSS(tx_subarrays, rx_subarrays_true, phi_current, theta_current, ...
                             phi_r, theta_r, lambda, d, N_sub);

    % 打印初始RSS值
    fprintf('\n重启点 %d, 初始方位角=%.2f, 初始极角=%.2f, 初始RSS: %.3f dB\n', restart_idx, phi_current, theta_current, rss_current);

    % 记录初始点的最佳值
    best_rss_local = rss_current;
    best_phi_local = phi_current;
    best_theta_local = theta_current;

    % ---- 新增: 初始化动量项 ----
    v_phi = 0;
    v_theta = 0;
    gamma = 0.9; % 动量参数
    % -------------------------

    % 执行梯度优化
    for iter = 1:max_iter
        % 记录当前值 [方位角, 极角]
        phi_history(iter) = phi_current;
        theta_history(iter) = theta_current;
        rss_history(iter) = rss_current;

        % --- 计算梯度 (dB) --- 
        phi_plus = phi_current + delta_angle;
        phi_minus = phi_current - delta_angle;
        phi_r_plus = phi_plus + 180; phi_r_plus = mod(phi_r_plus + 180, 360) - 180;
        phi_r_minus = phi_minus + 180; phi_r_minus = mod(phi_r_minus + 180, 360) - 180;
        theta_r_const = 180 - theta_current; theta_r_const = max(0, min(180, theta_r_const));
        % 使用 rx_subarrays_true
        rss_phi_plus = calculateRSS(tx_subarrays, rx_subarrays_true, phi_plus, theta_current, phi_r_plus, theta_r_const, lambda, d, N_sub);
        rss_phi_minus = calculateRSS(tx_subarrays, rx_subarrays_true, phi_minus, theta_current, phi_r_minus, theta_r_const, lambda, d, N_sub);
        gradient_phi = (rss_phi_plus - rss_phi_minus)/(2*delta_angle);

        theta_plus = theta_current + delta_angle;
        theta_minus = theta_current - delta_angle;
        theta_plus = max(0, min(180, theta_plus));
        theta_minus = max(0, min(180, theta_minus));
        phi_r_const = phi_current + 180; phi_r_const = mod(phi_r_const + 180, 360) - 180;
        theta_r_plus = 180 - theta_plus; theta_r_plus = max(0, min(180, theta_r_plus));
        theta_r_minus = 180 - theta_minus; theta_r_minus = max(0, min(180, theta_r_minus));
        % 使用 rx_subarrays_true
        rss_theta_plus = calculateRSS(tx_subarrays, rx_subarrays_true, phi_current, theta_plus, phi_r_const, theta_r_plus, lambda, d, N_sub);
        rss_theta_minus = calculateRSS(tx_subarrays, rx_subarrays_true, phi_current, theta_minus, phi_r_const, theta_r_minus, lambda, d, N_sub);
        if abs(theta_plus - theta_minus) > eps
            gradient_theta = (rss_theta_plus - rss_theta_minus)/(theta_plus - theta_minus);
        else
            gradient_theta = 0; 
        end
        
        if isinf(gradient_phi) || isnan(gradient_phi); gradient_phi = 0; fprintf('Inf/NaN grad_phi'); end
        if isinf(gradient_theta) || isnan(gradient_theta); gradient_theta = 0; fprintf('Inf/NaN grad_theta'); end
        
        fprintf('  Iter %d: Az Grad=%.2e, Pol Grad=%.2e, RSS=%.3f dB', iter, gradient_phi, gradient_theta, rss_current);

        % 梯度归一化
        grad_norm = sqrt(gradient_phi^2 + gradient_theta^2);
        if grad_norm > grad_zero_threshold
            norm_gradient_phi = gradient_phi / grad_norm;
            norm_gradient_theta = gradient_theta / grad_norm;
            grad_zero_count = 0; 
        else
            fprintf(', Grad < Threshold');
            grad_zero_count = grad_zero_count + 1;
            if iter >= min_iterations && grad_zero_count >= 5 
                fprintf('\nGradient near zero for %d iterations after %d min iterations, converging.\n', grad_zero_count, min_iterations);
                break;
            end
            fprintf('\n');
            continue; 
        end

        % --- 修改: 更新角度计算 (包含动量) ---
        % 计算原始更新量 (不含动量)
        delta_phi_update_raw = current_alpha_phi * norm_gradient_phi;
        delta_theta_update_raw = current_alpha_theta * norm_gradient_theta;

        % 计算包含动量的更新量
        v_phi = gamma * v_phi + delta_phi_update_raw; % 动量更新
        v_theta = gamma * v_theta + delta_theta_update_raw; % 动量更新

        % 尝试新位置
        phi_new = phi_current + v_phi; % 使用包含动量的更新
        theta_new = theta_current + v_theta; % 使用包含动量的更新
        % -----------------------------------
        phi_new = mod(phi_new + 180, 360) - 180;
        theta_new = max(0, min(180, theta_new));

        % 计算新接收端角度
        phi_r_new = phi_new + 180; phi_r_new = mod(phi_r_new + 180, 360) - 180;
        theta_r_new = 180 - theta_new; theta_r_new = max(0, min(180, theta_r_new));

        % 评估新位置 (使用 rx_subarrays_true)
        rss_new = calculateRSS(tx_subarrays, rx_subarrays_true, phi_new, theta_new, ...
                              phi_r_new, theta_r_new, lambda, d, N_sub);
        fprintf(', New RSS=%.3f dB', rss_new);

        % 更新策略
        if rss_new > rss_current + rss_improve_threshold 
            phi_current = phi_new;
            theta_current = theta_new;
            rss_current = rss_new;
            consecutive_decrease = 0; 
            fprintf(', Accepted.');
            current_alpha_phi = min(current_alpha_phi * 1.05, alpha_phi * 5); 
            current_alpha_theta = min(current_alpha_theta * 1.05, alpha_theta * 5);
        else
            consecutive_decrease = consecutive_decrease + 1;
            fprintf(', Rejected (Stuck %d)', consecutive_decrease);
            current_alpha_phi = max(current_alpha_phi * 0.7, alpha_phi / 10); 
            current_alpha_theta = max(current_alpha_theta * 0.7, alpha_theta / 10);
            
            if consecutive_decrease >= max_stuck_count
                fprintf('\nAttempting random perturbation... ');
                phi_current = phi_current + (rand-0.5) * perturb_angle;
                theta_current = theta_current + (rand-0.5) * perturb_angle;
                phi_current = mod(phi_current + 180, 360) - 180;
                theta_current = max(0, min(180, theta_current));
                
                phi_r = phi_current + 180; phi_r = mod(phi_r + 180, 360) - 180;
                theta_r = 180 - theta_current; theta_r = max(0, min(180, theta_r));
                % 使用 rx_subarrays_true
                rss_current = calculateRSS(tx_subarrays, rx_subarrays_true, phi_current, theta_current, phi_r, theta_r, lambda, d, N_sub);
                fprintf('New start Az=%.4f, Pol=%.4f, RSS=%.3f dB\n', phi_current, theta_current, rss_current);
                consecutive_decrease = 0; 
            end
        end
        fprintf('\n'); 

        % 记录本轮最佳位置
        if rss_current > best_rss_local
            best_rss_local = rss_current;
            best_phi_local = phi_current; 
            best_theta_local = theta_current;     
        end

        % 检查收敛条件 (基于角度变化和绝对RSS(dB)变化)
        if iter >= min_iterations
            prev_idx = max(iter-5, 1); 
            angle_change = sqrt(mean((phi_history(prev_idx:iter) - phi_history(iter)).^2) + ...
                                mean((theta_history(prev_idx:iter) - theta_history(iter)).^2)); 
            rss_abs_change = max(rss_history(prev_idx:iter)) - min(rss_history(prev_idx:iter));
            
            if angle_change < converge_threshold_angle && rss_abs_change < converge_threshold_rss_abs 
                fprintf('Convergence criteria met at iter %d (Angle=%.4f < %.4f, RSS_dB_change=%.2e < %.2e).\n', ...
                        iter, angle_change, converge_threshold_angle, rss_abs_change, converge_threshold_rss_abs);
                break; 
            end
        end
        
    end % End inner optimization loop

    % 更新全局最佳结果
    if best_rss_local > global_best_rss
        global_best_rss = best_rss_local;
        global_best_phi = best_phi_local; % Best Azimuth Overall
        global_best_theta = best_theta_local;     % Best Polar Overall
    end

    % 打印当前重启点的结果
    fprintf('Restart %d finished: Best Azimuth=%.4f°, Best Polar=%.4f°, Best RSS=%.3f dB\n', ...
        restart_idx, best_phi_local, best_theta_local, best_rss_local);
        
end % End restart loop

% 使用全局最佳结果作为最终结果
phi_opt = global_best_phi; % Final Azimuth
theta_opt = global_best_theta;     % Final Polar

fprintf('\nHigh-precision angle optimization result:\n');
fprintf('  Optimized Azimuth (phi): %.4f deg (True: %.4f deg, Error: %.4f deg)\n', phi_opt, phi_true, abs(mod(phi_opt - phi_true + 180, 360)-180));
fprintf('  Optimized Polar (theta): %.4f deg (True: %.4f deg, Error: %.4f deg)\n', theta_opt, theta_true, abs(theta_opt - theta_true));

%% 5. 结果可视化
figure('Position', [100, 100, 1200, 400]);

% 距离对比
subplot(1, 3, 1);
bar([D0, R0, R_opt]);
set(gca, 'XTickLabel', {'真实值', '感知输出', '对准估计'});
ylabel('距离 (m)');
title('距离估计对比');
grid on;
y_max = max([D0, R0, R_opt]);
y_min = min([D0, R0, R_opt]);
ylim([y_min*0.9, y_max*1.1]);

% 方位角对比
subplot(1, 3, 2);
bar([phi_true, phi0, phi_opt]);
set(gca, 'XTickLabel', {'真实值', '感知输出', '对准估计'});
ylabel('方位角 (phi, 度)');
title('方位角估计对比');
grid on;
y_max = max([phi_true, phi0, phi_opt]);
y_min = min([phi_true, phi0, phi_opt]);
ylim([y_min-2, y_max+2]); % Add margin

% 极角对比
subplot(1, 3, 3);
bar([theta_true, theta0, theta_opt]);
set(gca, 'XTickLabel', {'真实值', '感知输出', '对准估计'});
ylabel('极角 (theta, 度)');
title('极角估计对比');
grid on;
y_max = max([theta_true, theta0, theta_opt]);
y_min = min([theta_true, theta0, theta_opt]);
ylim([max(0, y_min-2), min(180, y_max+2)]); % Add margin within [0, 180]

% 优化过程可视化 (最后一个重启周期的)
figure('Position', [100, 550, 1200, 400]);
iter_actual = iter; % Actual iterations before break
if iter_actual == max_iter && grad_zero_count < 5 
   iter_actual = max_iter;
end
if iter_actual == 0; iter_actual = 1; end 
valid_indices = 1:iter_actual;

% 方位角优化过程
subplot(1, 3, 1);
plot(valid_indices, phi_history(valid_indices), 'o-', 'LineWidth', 1.5);
hold on;
plot([1, iter_actual], [phi_true, phi_true], 'r--', 'LineWidth', 1.5);
xlabel('迭代次数');
ylabel('方位角 (phi, 度)');
title('方位角优化过程 (Last Restart)');
legend('优化过程', '真实值', 'Location', 'best');
grid on;

% 极角优化过程
subplot(1, 3, 2);
plot(valid_indices, theta_history(valid_indices), 'o-', 'LineWidth', 1.5);
hold on;
plot([1, iter_actual], [theta_true, theta_true], 'r--', 'LineWidth', 1.5);
xlabel('迭代次数');
ylabel('极角 (theta, 度)');
title('极角优化过程 (Last Restart)');
legend('优化过程', '真实值', 'Location', 'best');
grid on;

% RSS优化过程 (dB)
subplot(1, 3, 3);
plot(valid_indices, rss_history(valid_indices), 'o-', 'LineWidth', 1.5); % Use plot for dB
xlabel('迭代次数');
ylabel('接收功率 (dB)');
title('信号强度优化过程 (Last Restart)');
grid on;

%% 结果汇总
fprintf('\n最终结果汇总:\n');
fprintf('参数          真实值      感知输出      对准估计      误差\n');
fprintf('----------------------------------------------------------------\n');
fprintf('距离 (m)      %.4f       %.4f         %.4f         %.4f\n', ...
    D0, R0, R_opt, abs(R_opt-D0));
fprintf('方位角 (phi, °) %.4f       %.4f         %.4f         %.4f\n', ...
    phi_true, phi0, phi_opt, abs(mod(phi_opt - phi_true + 180, 360)-180));
fprintf('极角 (theta, °) %.4f       %.4f         %.4f         %.4f\n', ...
    theta_true, theta0, theta_opt, abs(theta_opt - theta_true)); 

%% 6. 验证 RSS 在最终估计角度和真实角度的值
fprintf('\n--- RSS 验证 ---\n');

% 计算估计角度对应的接收端角度
phi_r_opt = phi_opt + 180;
phi_r_opt = mod(phi_r_opt + 180, 360) - 180;
theta_r_opt = 180 - theta_opt;
theta_r_opt = max(0, min(180, theta_r_opt));

% 计算真实角度对应的接收端角度
phi_r_true_calc = phi_true + 180;
phi_r_true_calc = mod(phi_r_true_calc + 180, 360) - 180;
theta_r_true_calc = 180 - theta_true;
theta_r_true_calc = max(0, min(180, theta_r_true_calc));

% 计算最终估计角度下的 RSS (确保 calculateRSS 设置与仿真时一致)
rss_at_estimated_angle = calculateRSS(tx_subarrays, rx_subarrays_true, phi_opt, theta_opt, ...
                                    phi_r_opt, theta_r_opt, lambda, d, N_sub);

% 计算真实角度下的 RSS (确保 calculateRSS 设置与仿真时一致)
rss_at_true_angle = calculateRSS(tx_subarrays, rx_subarrays_true, phi_true, theta_true, ...
                                 phi_r_true_calc, theta_r_true_calc, lambda, d, N_sub);

fprintf('RSS @ Estimated Angles (%.4f°, %.4f°): %.6f dB\n', phi_opt, theta_opt, rss_at_estimated_angle);
fprintf('RSS @ True Angles      (%.4f°, %.4f°): %.6f dB\n', phi_true, theta_true, rss_at_true_angle);

if rss_at_true_angle > rss_at_estimated_angle
    fprintf('注意：真实角度下的RSS高于估计角度下的RSS，优化器可能陷入局部最优。\n');
else
    fprintf('估计角度下的RSS不低于真实角度下的RSS。\n');
end 


%% 7. 诊断：绘制峰顶附近的 RSS 地形图
fprintf('\n--- 绘制 RSS 峰顶地形图 ---\n');

% 定义绘图范围和精度
plot_range_deg = 0.8; % 绘制范围 (+/- 真实角度)
plot_step_deg = 0.01; % 绘图步长

phi_range = (phi_true - plot_range_deg):plot_step_deg:(phi_true + plot_range_deg);
theta_range = (theta_true - plot_range_deg):plot_step_deg:(theta_true + plot_range_deg);

[PHI, THETA] = meshgrid(phi_range, theta_range);
RSS_map = zeros(size(PHI));

% 确保角度在有效范围内
THETA(THETA < 0) = 0;
THETA(THETA > 180) = 180;
PHI = mod(PHI + 180, 360) - 180;

% 计算网格上每个点的 RSS
fprintf('\nCalculating RSS map for diagnostics (may take a moment)...\n');
num_phi = length(phi_range);
num_theta = length(theta_range);

for i = 1:num_theta
    for j = 1:num_phi
        phi_t_map = PHI(i,j);
        theta_t_map = THETA(i,j);
        
        % 计算对应的接收端角度
        phi_r_map = phi_t_map + 180;
        phi_r_map = mod(phi_r_map + 180, 360) - 180;
        theta_r_map = 180 - theta_t_map;
        theta_r_map = max(0, min(180, theta_r_map));
        
        % 调用 calculateRSS (确保设置与仿真一致)
        RSS_map(i,j) = calculateRSS(tx_subarrays, rx_subarrays_true, phi_t_map, theta_t_map, ...
                                    phi_r_map, theta_r_map, lambda, d, N_sub);
        % 简单的进度指示
        if mod(j, round(num_phi/5)) == 0 && i==round(num_theta/2)
            fprintf('.'); 
        end
    end
end
fprintf('\nRSS map calculation complete.\n');

% 绘制 3D 表面图
figure('Name', 'RSS Peak Topography (3D)', 'Position', [200, 200, 800, 600]);
surf(PHI, THETA, RSS_map);
xlabel('方位角 (phi, 度)');
ylabel('极角 (theta, 度)');
zlabel('RSS (dB)');
title(sprintf('RSS 地形图 (中心: %.2f°, %.2f°)', phi_true, theta_true));
colorbar;
shading interp; % 平滑着色

% 寻找地图上的最大RSS值及其位置
[max_rss_map, idx_max_map] = max(RSS_map(:));
if ~isempty(idx_max_map) % 检查是否找到最大值
    [row_max, col_max] = ind2sub(size(RSS_map), idx_max_map(1)); % 取第一个最大值（可能多个点）
    phi_max_map = PHI(row_max, col_max);
    theta_max_map = THETA(row_max, col_max);
    fprintf('\n地图最大RSS: %.6f dB @ (%.4f°, %.4f°)\n', max_rss_map, phi_max_map, theta_max_map);
else
    fprintf('未能在地图上找到有效的最大RSS值。\n');
    phi_max_map = NaN;
    theta_max_map = NaN;
end

% 绘制 2D 等高线图
figure('Name', 'RSS Peak Topography (2D Contour)', 'Position', [1050, 200, 800, 600]);
contour(PHI, THETA, RSS_map, 50); % 绘制 50 条等高线
hold on;
plot(phi_true, theta_true, 'rx', 'MarkerSize', 12, 'LineWidth', 2); % 标记真实角度
plot(phi_opt, theta_opt, 'bo', 'MarkerSize', 10, 'LineWidth', 2);    % 标记优化结果
if ~isnan(phi_max_map)
    plot(phi_max_map, theta_max_map, 'g*', 'MarkerSize', 10, 'LineWidth', 2); % 标记地图峰值
    legend('RSS 等高线', '真实角度', '优化结果', '地图峰值', 'Location', 'best');
else
     legend('RSS 等高线', '真实角度', '优化结果', 'Location', 'best');
end
xlabel('方位角 (phi, 度)');
ylabel('极角 (theta, 度)');
title(sprintf('RSS 等高线图 (中心: %.2f°, %.2f°)', phi_true, theta_true));
colorbar;
grid on;
axis equal; % 保持横纵轴比例一致
