%% 扰动观察法（P&O）— 基于特性曲线表格实现
% 功能：先生成单二极管模型的 I-V/P-V 特性曲线，再在曲线上用变步长 P&O 搜索最大功率点
% 运行环境：MATLAB 或 Octave
% 说明：与 PandO_algorithm.m（逐点扰动版）互补，本脚本在预先算好的曲线上按索引搜索。
%       光伏模型参数统一来自 pv_params() + pv_current()。

%% 光伏组件参数（统一参数源）
pv = pv_params();
V_oc = pv.Voc_stc;
P_threshold = 2;   % 功率阈值(W)，用于变步长算法中判断是否远离功率峰值
V_threshold = 1;   % 电压阈值(V)，用于变步长算法中判断是否接近功率峰值

% 变步长系数定义
% 这些系数控制扰动观察法的步长调整策略，实现在不同功率变化区域采用不同步长以提高跟踪效率
k1 = 0.2;          % 远离峰值时的步长系数，采用较大步长快速接近最大功率点
k2 = 0.02;         % 接近峰值时的步长系数，采用较小步长精细调整避免振荡

%% 光伏特性曲线计算部分
% 通过统一单二极管模型计算光伏电池的I-V和P-V特性曲线，生成用于MPPT算法的基础数据
V_PV = linspace(0, V_oc, 1000).';  % 生成从0到开路电压的1000个均匀分布的电压点，转置为列向量
I_PV = max(pv_current(V_PV, pv), 0);  % 统一模型电流，保证非负
P = V_PV .* I_PV;   % 计算功率曲线，功率 = 电压 × 电流

% 将数据整理为向量形式，便于后续处理和分析（兼容 MATLAB 与 Octave）
V_vec = V_PV; I_vec = I_PV; P_vec = P;

%% MPPT算法初始化部分
% 设置算法迭代所需的初始参数和存储变量
N   = numel(V_vec);    % 获取电压点的总数
dV  = mean(diff(V_vec)); % 计算电压点的平均间隔，作为基础步长参考
iter = 1000;        % 设置最大迭代次数，确保算法能在有限步内收敛
idx = 1;            % 初始化电压索引，从第一个电压点开始搜索

% 预分配数组用于记录迭代过程数据，便于后续分析算法收敛性能
V_s = zeros(1,iter);  % 记录每次迭代的电压值
I_s = zeros(1,iter);  % 记录每次迭代的电流值
P_s = zeros(1,iter);  % 记录每次迭代的功率值

%% 变步长扰动观察法主循环
% 通过迭代调整工作点，跟踪最大功率点位置
for k = 1:iter
    % 获取当前工作点的电压和功率值
    V = V_vec(idx);   % 当前电压值
    P = P_vec(idx);   % 当前功率值
    % 获取相邻电压点的值用于计算变化量
    V_new = V_vec(idx + 1);  % 下一个电压点的值
    P_new = P_vec(idx + 1);  % 下一个电压点对应的功率值
    dP = P_new - P;    % 计算功率变化量，用于判断功率变化趋势
    dV = V_new - V;    % 计算电压变化量，即基础电压步长

    % 记录当前迭代的数据
    V_s(k) = V;    % 记录当前电压
    P_s(k) = P;    % 记录当前功率

    % 计算功率对电压的变化方向（斜率符号）
    % sign函数返回dP*dV的符号：正号表示同向变化，负号表示反向变化

    % 变步长策略：根据功率变化幅度和电压变化幅度动态调整步长
    if abs(dP) > P_threshold
        step = k1 * abs(dP);  % 当功率变化较大时，认为远离峰值，采用大步长快速接近
    elseif abs(dV) < V_threshold
        step = k2 * abs(dP);  % 当电压变化较小时，认为接近峰值，采用小步长精细调整
    else
        step = 0.1;  % 默认步长，用于中间状态
    end

    % 将电压步长转换为数据表格中的索引偏移量
    % 确保步长至少对应1个索引单位，避免除零错误
    step_idx = max(1, round(abs(step) / max(abs(dV), 1e-9)));
    % 根据功率斜率符号决定搜索方向：dP>0 则升压（右移），dP<0 则降压（左移）
    % 使用当前点 dP 的符号而非 sgn(dP*dV)，避免恒正导致的单向漂移
    step_idx = sign(P_new - P) * step_idx;
    % 计算新索引位置
    idx_new = idx + step_idx;
    % 限制索引在有效范围内[1, N-1]，避免越界
    idx = min(max(idx_new, 1), N-1);
end

%% 输出最终结果
% 算法收敛后，获取最大功率点对应的电压和功率值
V_mppt = V_vec(idx);  % 最大功率点电压
P_mppt = P_vec(idx);  % 最大功率点功率

%% 结果可视化部分
% 绘制I-V特性曲线，展示光伏电池的电流电压关系
figure; 
plot(V_vec, I_vec, '-b');  % 绘制I-V曲线，蓝色实线
xlabel('Voltage (V)');  % X轴标签：电压
ylabel('Current (A)');  % Y轴标签：电流
title('I–V');          % 图表标题：I-V特性曲线
hold on; 
plot(V_mppt, I_vec(idx), 'ro');  % 在曲线上标记最大功率点，红色圆圈

% 绘制P-V特性曲线，展示光伏电池的功率电压关系
figure; 
plot(V_vec, P_vec, '-r');  % 绘制P-V曲线，红色实线
xlabel('Voltage (V)');  % X轴标签：电压
ylabel('Power (W)');    % Y轴标签：功率
title('P–V');          % 图表标题：P-V特性曲线
hold on; 
plot(V_mppt, P_mppt, 'ko');  % 标记最大功率点，黑色圆圈

% 绘制功率迭代过程，展示算法收敛性能
figure; 
plot(P_s,'LineWidth',1.2);  % 绘制功率随迭代次数的变化曲线
xlabel('迭代过程');          % X轴标签：迭代次数
ylabel('P (W)');            % Y轴标签：功率

% 绘制电压迭代过程，展示工作点电压的变化
figure; 
plot(V_s,'LineWidth',1.2);  % 绘制电压随迭代次数的变化曲线  
xlabel('迭代过程');          % X轴标签：迭代次数
ylabel('V_{pv} (V)');       % Y轴标签：光伏电压
