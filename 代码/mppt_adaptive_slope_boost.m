%% pv_mppt_boost_slopePO.m
% 单文件：PV 单二极管模型 + 理想 Boost 平均模型 + 变步长(基于 dP/dV) P&O MPPT
% 运行环境：MATLAB/Octave (fzero在Octave需optim包)
clear; clc; close all;  % 清除工作区变量、命令窗口和所有图形窗口

% ==== PV 组件与物理参数（可按实际组件微调）====
pv.Ns        = 60;      % 串联片数，光伏组件中单体电池的串联数量
pv.Rs        = 0.35;    % 串联电阻 [Ohm]，模拟电池内部和连接线路的电阻损耗
pv.Rsh       = 500;     % 并联电阻 [Ohm]，模拟电池的漏电流路径
pv.n         = 1.3;     % 理想因子，反映二极管非理想特性，典型值1~2
pv.Isc_stc   = 8.7;     % STC短路电流 [A]，标准测试条件下的短路电流
pv.Voc_stc   = 37.0;    % STC开路电压 [V]，标准测试条件下的开路电压
pv.alpha_Isc = 0.005;   % 电流温度系数 [A/°C]，短路电流随温度的变化率
pv.Eg        = 1.12;    % 硅材料禁带宽度 [eV]，影响反向饱和电流的温度特性

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

% ===== 局部函数定义 =====

function [Vsol, Isol] = solve_pv_operating_point(D, Rload, G, T, params)
% 功能：求解给定占空比下的光伏工作点电压和电流
% 输入：D-占空比, Rload-负载电阻, G-辐照度, T-温度, params-光伏参数
% 输出：Vsol-工作点电压, Isol-工作点电流

% 计算Boost变换器在光伏侧看到的等效输入电阻
% 理想Boost变换器关系：Rin = Rload * (1-D)^2
Rin = Rload * (1 - D)^2;
if Rin <= 0, Rin = 1e-6; end  % 防零除保护

% 定义残差方程：光伏输出电流等于负载电流
% f(V) = Ipv(V) - V/Rin = 0
f = @(V) pv_current_given_voltage(max(V,0), G, T, params) - max(V,0)/Rin;

% 估算开路电压作为搜索区间上限（简化温度修正）
Voc_est = params.Voc_stc * (1 + 0.0*(T-25)/25); % 简化处理，实际应加温度系数

% TODO: 开路电压温度修正
% 问题：未考虑温度对开路电压的影响，实际Voc随温度升高而降低
% 建议：添加温度系数，如 Voc_est = params.Voc_stc + beta*(T-25)，beta ≈ -0.003/°C

V0 = 0.7*Voc_est;  % 初始猜测值，设为估计开路电压的70%

% TODO: fzero初始猜测优化
% 问题：固定比例初始值在极端条件下（如低辐照度）可能不理想
% 建议：基于当前工作条件动态调整初始猜测值

% 使用fzero求解非线性方程
try
    % 优先使用区间搜索，提高收敛可靠性
    Vsol = fzero(f, [0, max(1e-6, 1.2*Voc_est)]);
catch
    % 区间搜索失败时使用单点搜索
    Vsol = fzero(f, V0);
end
Vsol = max(Vsol, 0);  % 确保电压非负
Isol = pv_current_given_voltage(Vsol, G, T, params);  % 计算对应电流
end

function I = pv_current_given_voltage(V, G, T, params)
% 功能：基于单二极管模型计算给定电压下的光伏输出电流
% 输入：V-电压, G-辐照度, T-温度, params-光伏参数
% 输出：I-电流

% 物理常数定义
q = 1.602176634e-19;  % 元电荷 [C]
k = 1.380649e-23;     % 玻尔兹曼常数 [J/K]
T_K = T + 273.15;     % 转换为绝对温度 [K]
Vt  = params.Ns * k*T_K / q;  % 热电压，考虑串联电池数
n   = params.n;               % 理想因子

% 计算光生电流 Iph（考虑辐照度和温度影响）
Iph_stc = params.Isc_stc;  % STC条件下的光生电流
Iph = Iph_stc*(G/1000) + params.alpha_Isc*(T - 25);  % 实际条件下的光生电流

% 计算反向饱和电流 I0（考虑温度影响）
Vt_stc = params.Ns * k*(25+273.15)/q;  % STC条件下的热电压
% 由STC条件反推I0_stc：Isc ≈ Iph ≈ I0*exp(Voc/(nVt))
I0_stc = Iph_stc / (exp(params.Voc_stc/(n*Vt_stc)) - 1 + 1e-12);

% TODO: 数值稳定性改进 - 指数溢出保护
% 问题：当指数参数过大时，exp()函数可能产生Inf，特别是在高电压条件下
% 建议：添加指数参数限制，如 exp(min(x, 50))

Eg = params.Eg;  % 材料禁带宽度

% I0的温度修正公式，考虑温度对载流子浓度和扩散系数的影响
I0 = I0_stc * (T_K/(25+273.15)).^3 .* ...
    exp( (q*Eg/k) * (1/(25+273.15) - 1/T_K) / (n*params.Ns) );

% 提取串联和并联电阻参数
Rs  = params.Rs;   % 串联电阻
Rsh = params.Rsh;  % 并联电阻

% 定义单二极管模型的隐式方程：f(I) = 0
fI = @(I) Iph - I0.*(exp((V + I.*Rs)./(n*Vt)) - 1) - (V + I.*Rs)./Rsh - I;

% TODO: 电流计算数值稳定性增强
% 问题：指数函数在高压条件下可能数值溢出，导致求解失败
% 建议：在指数函数中添加参数限制，提高数值稳定性

% 设置合理的初始猜测值
I_init = max(0, min(Iph, Iph - V/max(Rsh,1e-6)));

% TODO: fzero异常处理增强
% 问题：当前的try-catch可能无法处理所有异常情况，特别是数值不稳定的情况
% 建议：添加更详细的错误处理，包括多个备用求解策略

% 求解非线性方程得到电流值
try
    I = fzero(fI, I_init);
catch
    % 单点搜索失败时使用区间搜索
    I = fzero(fI, [0, max(0.1, Iph*1.5)]);
end
I = max(I, 0);  % 确保电流非负
end

function [Vvec, Ivec] = sweep_pv_iv(G, T, params)
% 功能：扫描生成完整I-V特性曲线
% 输入：G-辐照度, T-温度, params-光伏参数
% 输出：Vvec-电压向量, Ivec-电流向量

Vmax = params.Voc_stc*1.15;  % 设置扫描电压上限（略大于开路电压）
Vvec = linspace(0, Vmax, 200);  % 生成200个电压点
Ivec = zeros(size(Vvec));       % 预分配电流数组

% 遍历所有电压点计算对应电流
for i = 1:numel(Vvec)
    Ivec(i) = pv_current_given_voltage(Vvec(i), G, T, params);
end
end
