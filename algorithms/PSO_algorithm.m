% 粒子群算法
% 参数设置
num_particles = 30;   % 粒子数量
num_dimensions = 3;   % 搜索空间维度
max_iter = 100;       % 最大迭代次数
w = 0.5;              % 惯性权重
c1 = 1.5;             % 个体加速常数
c2 = 2.0;             % 社会加速常数
lb = -10;             % 搜索空间下界
ub = 10;              % 搜索空间上界
tolerance = 1e-5

% 初始化粒子群
particles = rand(num_particles, num_dimensions) * (ub - lb) + lb;
velocities = zeros(num_particles, num_dimensions);
p_best = particles; % 个体最优位置
p_best_scores = inf(num_particles, 1); % 个体最优得分
g_best = particles(1, :); % 全局最优位置
g_best_score = inf;       % 全局最优得分

% PSO主循环
for iter = 1:max_iter
    % 评估每个粒子
    for i = 1:num_particles
        % 计算粒子适应度
        current_score = evaluate(particles(i, :));
        
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
