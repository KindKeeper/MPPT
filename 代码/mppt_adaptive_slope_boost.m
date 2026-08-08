%% pv_mppt_boost_slopePO.m
% 单文件：PV 单二极管模型 + 理想 Boost 平均模型 + 变步长(基于 dP/dV) P&O MPPT
% 运行环境：MATLAB 或 Octave（fzero 为两者内置，无需额外包）
% 说明：依赖的三个辅助函数已拆分为独立文件，位于同目录：
%       solve_pv_operating_point.m / pv_current_given_voltage.m / sweep_pv_iv.m
%       （拆分后 MATLAB 与 Octave 均可运行，函数也便于复用）
%       光伏模型参数统一来自 pv_params() + pv_current()。
clear; clc; close all;  % 清除工作区变量、命令窗口和所有图形窗口

% ==== PV 组件与物理参数（统一参数源，可按实际组件微调 pv_params.m）====
pv = pv_params();

% 环境参数定义
G  = 800;               % 辐照度 [W/m^2]，当前环境的光照强度
Tc = 35;                % 组件温度 [°C]，光伏板的工作温度

% ==== 负载与Boost平均模型 ====
Rload = 20;             % 直流负载电阻 [Ohm]，Boost变换器输出侧的等效负载
Dmin = 0.02; Dmax = 0.95;   % 占空比安全边界，避免极端工况导致系统不稳定
D    = 0.5;                 % 初始占空比，MPPT算法的起始工作点

% 说明：在理想连续导通(CCM)平均模型下，源侧看到的等效阻抗
% Rin = Rload * (1 - D)^2
% 工作点满足 Ipv(V) = V / Rin。给定D可解得 Vpv。

% ==== MPPT（基于斜率 dP/dV 的自适应扰动）参数 ====
Nstep      = 300;      % 最大迭代步数，控制算法运行时间
Kgain      = 1e-3;     % 步长增益系数，将功率斜率转换为步长的比例因子
step_min   = 2e-4;     % 最小步长，在最大功率点附近使用，减少振荡
step_max   = 1.5e-2;   % 最大步长，在远离最大功率点时使用，加快收敛
slope_eps  = 1e-3;     % 斜率死区阈值，小于此值认为已接近最大功率点
slope_LP   = 0.7;      % 斜率低通滤波系数，平滑斜率估计，减少噪声影响

% 数据记录数组预分配，用于存储迭代过程数据
Dlog = zeros(1,Nstep);  % 占空比记录
Vlog = zeros(1,Nstep);  % 光伏电压记录
Ilog = zeros(1,Nstep);  % 光伏电流记录
Plog = zeros(1,Nstep);  % 光伏功率记录

% 为首步初始化一个参考点
[Vk, Ik] = solve_pv_operating_point(D, Rload, G, Tc, pv);  % 求解初始工作点
Pk = Vk * Ik;              % 计算初始功率
dPdV_est_f = 0;           % 初始化滤波后的功率斜率估计值

% TODO: 低通滤波器初始化优化
% 问题：滤波器初始化为0可能导致初始收敛慢，第一次迭代的滤波值不准确
% 建议：在循环前先进行一次斜率测量，用实际测量值初始化滤波器

%% MPPT主循环 - 基于功率斜率的自适应扰动观察法
for k = 1:Nstep
    % 1) 记录当前迭代数据
    Dlog(k) = D;    % 记录当前占空比
    Vlog(k) = Vk;   % 记录当前电压
    Ilog(k) = Ik;   % 记录当前电流
    Plog(k) = Pk;   % 记录当前功率

    % 2) 预测下一个操作点以估计斜率 dP/dV（数值差分法）
    D_probe = min(max(D + 1e-3, Dmin), Dmax);     % 施加微小扰动用于斜率估计

    % TODO: 斜率估计步长优化
    % 问题：固定扰动步长1e-3在不同工作条件下可能不是最优选择
    % 建议：基于当前电压动态调整扰动步长，如使用百分比方式

    [V_probe, I_probe] = solve_pv_operating_point(D_probe, Rload, G, Tc, pv);
    P_probe = V_probe * I_probe;  % 计算扰动后的功率

    % 计算电压和功率的变化量
    dV   = (V_probe - Vk);
    dP   = (P_probe - Pk);
    % 计算功率对电压的导数（斜率），添加防除零保护
    dPdV = (abs(dV) > 1e-9) * (dP / (dV + (abs(dV)<=1e-9))) ;
    % 使用一阶低通滤波器平滑斜率估计，减少测量噪声影响
    dPdV_est_f = slope_LP*dPdV_est_f + (1 - slope_LP)*dPdV;

    % 3) 依据 dP/dV 确定搜索方向与自适应步长
    if abs(dPdV_est_f) < slope_eps
        % 接近最大功率点区域，斜率很小
        dir = 0;                    % 方向为零，减少不必要的扰动
        step = step_min;            % 使用最小步长，降低功率振荡
    else
        % 远离最大功率点，根据斜率方向确定搜索方向
        dir = sign(dPdV_est_f);     % 斜率正负决定电压调整方向
        % 自适应步长：斜率越大步长越大，但限制在[min, max]范围内
        step = min(step_max, max(step_min, Kgain*abs(dPdV_est_f)));
    end

    % 4) 更新占空比 (Boost变换器特性：D↑ -> Vpv↓；D↓ -> Vpv↑)
    % 注意：希望V增加时(dPdV>0)，应减小D；希望V减小时(dPdV<0)，应增大D
    D = D - dir*step;              % 根据方向调整占空比
    D = min(max(D, Dmin), Dmax);   % 限制占空比在安全范围内

    % 5) 用新占空比求解新的工作点
    [Vk, Ik] = solve_pv_operating_point(D, Rload, G, Tc, pv);
    Pk = Vk * Ik;  % 更新功率值
end

%% ==== 结果可视化部分 ====
% 生成当前环境条件下的完整I-V特性曲线
[Vvec, Ivec] = sweep_pv_iv(G, Tc, pv);
Pvec = Vvec .* Ivec;  % 计算P-V曲线

% 绘制I-V特性曲线及MPPT最终工作点
figure;
plot(Vvec, Ivec, 'LineWidth',1.5); grid on; hold on;
plot(Vlog(end), Ilog(end), 'o', 'MarkerSize',6);
xlabel('V_{pv} (V)'); ylabel('I_{pv} (A)');
title('PV I–V (current G,T) & MPPT result');
legend('I–V','MPPT point','Location','best');

% 绘制P-V特性曲线及MPPT最终工作点
figure;
plot(Vvec, Pvec, 'LineWidth',1.5); grid on; hold on;
plot(Vlog(end), Plog(end), 'o', 'MarkerSize',6);
xlabel('V_{pv} (V)'); ylabel('P_{pv} (W)');
title('PV P–V (current G,T) & MPPT result');
legend('P–V','MPPT point','Location','best');

% 绘制功率收敛过程
figure;
plot(Plog,'LineWidth',1.2); grid on;
xlabel('Iteration'); ylabel('P (W)');
title('Convergence of Power');

% 绘制电压变化过程
figure;
plot(Vlog,'LineWidth',1.2); grid on;
xlabel('Iteration'); ylabel('V_{pv} (V)');
title('PV Voltage vs. Iteration');

% 绘制占空比变化过程
figure;
plot(Dlog,'LineWidth',1.2); grid on;
xlabel('Iteration'); ylabel('D');
title('Duty Ratio vs. Iteration');

