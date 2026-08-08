%% 粒子群算法（PSO）搜索MPPT最优电压
% 功能：使用粒子群算法搜索光伏最大功率点对应的最优电压
% 适应度函数：目标为使功率 P = V*I 最大（PSO默认求最小，故取负）
% 运行环境：MATLAB 或 Octave
% 说明：本脚本为独立完整示例，直接运行即可看到搜索结果。
%       光伏模型参数统一来自 pv_params() + pv_current()。

clear; clc; close all;

% ===== 光伏组件参数（统一参数源） =====
pv = pv_params();
V_oc = pv.Voc_stc;

% 单二极管模型电流-电压关系（统一模型入口）
pv_iv = @(V) pv_current(V, pv);

% 适应度函数：功率取负（PSO求最小），元素对应点乘
fitness = @(V) -(V .* pv_iv(V));

% ===== PSO 参数设置 =====
num_particles = 30;   % 粒子数量
num_dimensions = 1;   % 搜索空间维度（最优电压）
max_iter = 100;       % 最大迭代次数
w = 0.5;              % 惯性权重
c1 = 1.5;             % 个体加速常数
c2 = 2.0;             % 社会加速常数
lb = 0;               % 搜索空间下界
ub = V_oc;            % 搜索空间上界
tolerance = 1e-5;     % 收敛容差

% ===== 初始化粒子群 =====
particles = rand(num_particles, num_dimensions) * (ub - lb) + lb;
velocities = zeros(num_particles, num_dimensions);
p_best = particles; % 个体最优位置
p_best_scores = inf(num_particles, 1); % 个体最优得分
g_best = particles(1, :); % 全局最优位置
g_best_score = inf;       % 全局最优得分

% ===== PSO主循环 =====
for iter = 1:max_iter
    % 评估每个粒子
    for i = 1:num_particles
        % 计算粒子适应度
        current_score = fitness(particles(i, :));
        
        % 更新个体最优
        if current_score < p_best_scores(i)
            p_best_scores(i) = current_score;
            p_best(i, :) = particles(i, :);
        end
        
        % 更新全局最优
        if current_score < g_best_score
            g_best_score = current_score;
            g_best = particles(i, :);
        end
    end
 
    % 更新粒子速度和位置
    for i = 1:num_particles
        % 标准速度更新公式
        velocities(i, :) = w * velocities(i, :) ...
            + c1 * rand * (p_best(i, :) - particles(i, :)) ...
            + c2 * rand * (g_best - particles(i, :));
        
        % 速度边界限制
        velocities(i, :) = max(min(velocities(i, :), ub), lb);
        
        % 位置更新
        particles(i, :) = particles(i, :) + velocities(i, :);
        
        % 位置边界限制
        particles(i, :) = max(min(particles(i, :), ub), lb);
    end
    if std(p_best_scores) < tolerance
       break; % 群体收敛时提前终止
    end
end

% ===== 结果展示 =====
V_mppt = g_best(1);
I_mppt = pv_iv(V_mppt);
P_mppt = V_mppt * I_mppt;
fprintf('PSO找到的最优电压: %.2f V\n', V_mppt);
fprintf('对应电流: %.2f A, 最大功率: %.2f W\n', I_mppt, P_mppt);
