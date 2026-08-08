%% 光伏MPPT算法性能对比：扰动观察法(P&O) vs 粒子群算法(PSO)
% 功能：在统一单二极管光伏模型下，分别运行P&O和PSO寻找最大功率点，
%       对比两者收敛结果与迭代开销，并输出对比结论。
% 运行环境：MATLAB 或 Octave
% 说明：本脚本为独立完整的对比演示，便于学习和验证。
%       光伏模型参数统一来自 pv_params() + pv_current()。

clear; clc; close all;

%% ================== 光伏模型（统一单二极管模型）==================
pv = pv_params();
V_oc = pv.Voc_stc;

% 单二极管模型电流-电压关系（统一模型入口）
pv_equation = @(V) pv_current(V, pv);


%% ================= 方法一：扰动观察法(P&O) =================
% 通过扰动电压，观察功率变化方向逐步逼近最大功率点
delta = 0.05;    % 扰动步长 [V]
V = 1;           % 初始电压 [V]
P_po_hist = [];
for iter = 1:5000
    I = pv_equation(V);
    P = V * I;
    P_po_hist(end+1) = P;
    % 尝试正向扰动
    V_new = V + delta;
    I_new = pv_equation(V_new);
    P_new = V_new * I_new;
    if P_new > P
        V = V_new;   % 功率上升，保持方向
    else
        V = V - delta;  % 功率下降，反向扰动
    end
    V = min(max(V, 0), V_oc);   % 限幅
    if abs(P_new - P) < 1e-3
        break;
    end
end
V_po = V; I_po = pv_equation(V_po); P_po = V_po * I_po;

%% ================= 方法二：粒子群算法(PSO) =================
% 粒子群优化搜索电压，使功率最大（适应度 = 功率）
num_particles = 30;
max_iter = 50;
w = 0.5;  c1 = 1.5;  c2 = 2.0;
lb = 0;  ub = 36.6;   % 搜索电压范围
particles = lb + (ub-lb).*rand(num_particles,1);
velocities = zeros(num_particles,1);
p_best = particles;
p_best_score = arrayfun(@(v) v*pv_equation(v), particles);
[g_best_score, gi] = max(p_best_score);
g_best = particles(gi);
for iter = 1:max_iter
    for i = 1:num_particles
        score = particles(i) * pv_equation(particles(i));
        if score > p_best_score(i)
            p_best_score(i) = score;
            p_best(i) = particles(i);
        end
        if score > g_best_score
            g_best_score = score;
            g_best = particles(i);
        end
    end
    for i = 1:num_particles
        r1 = rand(); r2 = rand();
        velocities(i) = w*velocities(i) + c1*r1*(p_best(i)-particles(i)) + c2*r2*(g_best-particles(i));
        particles(i) = max(lb, min(ub, particles(i) + velocities(i)));
    end
end
V_pso = g_best;  I_pso = pv_equation(V_pso);  P_pso = V_pso * I_pso;

%% ================= 结果展示与对比 =================
% 绘制P-V曲线并标出两种算法的收敛点
V_sweep = linspace(0, V_oc, 500);
P_sweep = V_sweep .* pv_equation(V_sweep);
figure;
plot(V_sweep, P_sweep, '-', 'LineWidth', 1.5); hold on;
plot(V_po, P_po, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
plot(V_pso, P_pso, 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('V (V)'); ylabel('P (W)'); grid on;
legend('P-V曲线','P&O 结果','PSO 结果','Location','best');
title('MPPT 算法对比');

fprintf('—— 对比结果 ——\n');
fprintf('P&O : V=%.2fV, I=%.2fA, P=%.2fW\n', V_po, I_po, P_po);
fprintf('PSO : V=%.2fV, I=%.2fA, P=%.2fW\n', V_pso, I_pso, P_pso);
fprintf('P&O 迭代次数: %d\n', numel(P_po_hist));
fprintf('PSO 迭代代数: %d\n', max_iter);