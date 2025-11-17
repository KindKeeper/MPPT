%% 光伏MPPT算法实现与性能测试
% 本文件包含两种MPPT算法实现：扰动观察法(P&O)和粒子群算法(PSO)
% 以及代码质量测试方案
% 功能：比较不同MPPT算法的性能，并提供代码质量评估工具

%% TODO: 拆分方案建议
% 方案1：拆分为三个独立文件
%   P_O_algorithm.m - 扰动观察法实现
%   PSO_algorithm.m - 粒子群算法实现(包含evaluate函数)
%   test_code_quality.m - 代码质量测试脚本

% 方案2：拆分为四个文件
%   P_O_algorithm.m - 扰动观察法
%   PSO_main.m - 粒子群主程序
%   evaluate.m - 粒子群评估函数
%   test_code_quality.m - 测试脚本

%% ================== 扰动观察法(P&O)实现 ==================
% 功能：通过扰动观察法寻找光伏电池的最大功率点
% 原理：通过小步长扰动电压，观察功率变化方向，逐步逼近最大功率点

% TODO: 优化建议1 - 添加步长自适应机制
%   可考虑: delta = delta0 * abs(dP/dV); % 根据功率变化率调整步长
% TODO: 优化建议2 - 添加环境突变检测
%   可考虑: if abs(I_L - prev_I_L) > threshold
%              reset_algorithm(); % 环境突变时重置算法
%          end

% 光伏电池参数定义
I_L = 10;     % 光生电流 [A]，光伏电池在光照下产生的电流
I_o = 1e-7;   % 反向饱和电流 [A]，反映PN结反向特性的参数
q = 1.6e-19;  % 元电荷 [C]，基本电荷常数1.6×10^-19库仑
T = 298;      % 室温 [K]，工作温度298K=25°C
n = 1.5;      % 理想因子，二极管非理想程度参数，典型值1-2
k = 1.38e-23; % 玻尔兹曼常数 [J/K]，热力学统计常数
delta = 0.01; % 扰动步长 [V]，电压扰动量，影响收敛速度和精度
V = 10;       % 初始电压 [V]，算法开始搜索的起始电压值

% P&O算法主循环 - 通过迭代扰动电压寻找最大功率点
for iter = 1:1000 
    % 计算当前工作点的电流和功率
    % 使用单二极管模型计算输出电流：I = I_L - I_o*(exp(qV/nkT) - 1)
    I = I_L - I_o*(exp((q*V)/(n*k*T)) - 1);
    % 计算当前功率：P = V × I
    P = V * I;
    
    % 电压扰动：在当前电压基础上施加小扰动
    V_new = V + delta;  % 新电压点
    % 计算新电压点对应的电流值
    I_new = I_L - I_o*(exp((q*V)/(n*k*T)) - 1);
    % 计算新电压点对应的功率值
    P_new = V_new * I_new;
    
    % 判断功率变化方向，决定下一步扰动方向
    if P_new > P
        % 功率增加，保持当前扰动方向（继续正向扰动）
        delta = delta; % 继续正向扰动
    else
        % 功率减小，改变扰动方向（反向扰动）
        delta = -delta; % 反向扰动
    end
    
    % 更新电压：应用扰动得到新的工作电压
    V = V + delta;
    
    % 终止条件: 功率变化小于阈值，认为已收敛到最大功率点
    if abs(P_new - P) < 1e-3
        break;  % 功率变化足够小，退出循环
    end
end

% 输出MPPT结果：显示找到的最大功率点电压和功率
disp(['MPP电压: ', num2str(V), 'V, 功率: ', num2str(P), 'W']);

%% ================== 粒子群算法(PSO)实现 ==================
% 功能：使用粒子群优化算法寻找多维空间中的最优解
% 原理：模拟鸟群觅食行为，通过个体和群体经验指导搜索方向

% TODO: 优化建议1 - 添加惯性权重自适应
%   可考虑: w = w_max - (w_max-w_min)*iter/max_iter;
% TODO: 优化建议2 - 添加早停机制
%   可考虑: if std(p_best_scores) < tolerance
%              break; % 群体收敛时提前终止
%           end

% PSO参数设置
num_particles = 30;   % 粒子数量，影响搜索的全面性和计算复杂度
num_dimensions = 5;   % 搜索空间维度，问题的变量个数
max_iter = 100;       % 最大迭代次数，控制算法运行时间
w = 0.5;              % 惯性权重，平衡全局和局部搜索能力（0.4-0.9）
c1 = 1.5;             % 个体加速常数，控制个体经验的影响力
c2 = 2.0;             % 社会加速常数，控制群体经验的影响力
lb = -10;             % 搜索空间下界，变量最小值约束
ub = 10;              % 搜索空间上界，变量最大值约束

% 初始化粒子群：在搜索空间内随机生成初始粒子位置
particles = rand(num_particles, num_dimensions) * (ub - lb) + lb;
% 初始化粒子速度：开始时速度为零
velocities = zeros(num_particles, num_dimensions);
% 初始化个体最优位置：开始时每个粒子的位置就是其个体最优
p_best = particles; % 个体最优位置矩阵
% 初始化个体最优得分：设为无穷大，便于后续比较
p_best_scores = inf(num_particles, 1); % 个体最优得分向量
% 初始化全局最优位置：暂时设为第一个粒子的位置
g_best = particles(1, :); % 全局最优位置向量
% 初始化全局最优得分：设为无穷大
g_best_score = inf;       % 全局最优得分

% PSO主循环：通过迭代更新粒子位置寻找最优解
for iter = 1:max_iter
    % 评估每个粒子的适应度
    for i = 1:num_particles
        % 计算粒子适应度（当前解的优劣程度）
        current_score = evaluate(particles(i, :));
        
        % 更新个体最优：如果当前解优于个体历史最优
        if current_score < p_best_scores(i)
            p_best_scores(i) = current_score;  % 更新个体最优得分
            p_best(i, :) = particles(i, :);   % 更新个体最优位置
        end
        
        % 更新全局最优：如果当前解优于全局历史最优
        if current_score < g_best_score
            g_best_score = current_score;  % 更新全局最优得分
            g_best = particles(i, :);      % 更新全局最优位置
        end
    end
 
    % 更新粒子速度和位置：根据个体和群体经验调整搜索方向
    for i = 1:num_particles
        % 速度更新公式：新速度 = 惯性部分 + 个体认知部分 + 社会认知部分
        velocities(i, :) = w * velocities(i, :) ...          % 惯性部分
            + c1 * rand * (p_best(i, :) - particles(i, :)) ... % 个体认知
            + c2 * rand * (g_best - particles(i, :));         % 社会认知
        
        % 速度边界限制：防止速度过大导致搜索不稳定
        velocities(i, :) = max(min(velocities(i, :), ub), lb);
        
        % 位置更新：新位置 = 原位置 + 速度
        particles(i, :) = particles(i, :) + velocities(i, :);
        
        % 位置边界限制：确保粒子在搜索空间内
        particles(i, :) = max(min(particles(i, :), ub), lb);
    end
end

%% TODO: 以下代码可能冗余，建议检查后移除
% 更新粒子速度（此段代码可能重复，建议检查）
velocities(i, :) = w * velocities(i, :) ...
    + c1 * rand * (p_best(i, :) - particles(i, :)) ...
    + c2 * rand * (g_best - particles(i, :));

% 适应度评估函数
function score = evaluate(position)
    % 功能：评估粒子位置的适应度（解的质量）
    % 输入：position - 粒子位置向量
    % 输出：score - 适应度得分（越小越好）
    
    % TODO: 优化建议 - 添加实际问题适配
    % 当前为测试函数，实际应用中需替换为光伏系统模型
    % 这里使用简单的球函数作为测试：f(x) = -sum(x^2)
    % 实际光伏应用中应改为：f(D) = -P(D)，即寻找使功率最大的占空比
    score = -sum(position.^2); 
end

%% ================== 代码质量测试方案 ==================
% 功能：提供算法性能测试和代码质量分析工具
% 包括运行时间测试、CPU时间测试、静态代码分析和性能剖析

% TODO: 优化建议 - 封装为独立测试函数
%   可考虑: function run_tests()
%              % 测试代码放在这里
%           end

% 运行时间测试：测量算法实际执行时间（墙钟时间）
timerVal = tic;  % 开始计时
% 执行算法 (实际测试时需取消注释)
% P_O_algorithm;
% PSO_algorithm;
elapsedTime = toc(timerVal);  % 结束计时并计算耗时
disp(['总运行时间: ', num2str(elapsedTime), '秒']);

% CPU时间测试：测量算法消耗的CPU时间（更准确的性能指标）
t_Start = cputime;  % 记录开始CPU时间
% 执行算法 (实际测试时需取消注释)
% P_O_algorithm;
% PSO_algorithm;
t_End = cputime - t_Start;  % 计算CPU时间差
disp(['CPU时间: ', num2str(t_End), '秒']);

% 代码静态分析：检查代码潜在问题和改进建议
% TODO: 优化建议 - 添加自动化分析报告生成
checkcode('当前文件名.m');  % 执行MATLAB代码检查器
info = checkcode('当前文件名.m', '-id');  % 获取详细的检查信息
disp([info.message]);  % 显示检查结果信息

% 性能剖析：分析代码中各部分的执行时间分布
profile on  % 开启性能剖析
% 执行算法 (实际测试时需取消注释)
% P_O_algorithm;
% PSO_algorithm;
p = profile('info');  % 获取剖析信息
save myprofiledata p;  % 保存剖析数据供后续分析
profile viewer  % 打开性能剖析查看器
