"""
pso_algorithm.py — 粒子群算法（PSO）搜索 MPPT 最优电压
对应 PSO_algorithm.m：粒子群优化电压使功率最大（适应度=功率，PSO 求最小故取负）。
模型统一来自 pv_model.pv_params() + pv_current()。
运行：uv run python pso_algorithm.py
"""
import numpy as np
from pv_model import pv_params, pv_current

rng = np.random.default_rng(42)   # 固定随机种子，保证结果可复现

pv = pv_params()
V_oc = pv["Voc_stc"]
pv_iv = lambda V: pv_current(V, pv)

fitness = lambda V: -(V * pv_iv(V))   # PSO 求最小，故功率取负

# PSO 参数
num_particles = 30
num_dimensions = 1
max_iter = 100
w, c1, c2 = 0.5, 1.5, 2.0
lb, ub = 0, V_oc
tolerance = 1e-5

particles = rng.uniform(lb, ub, (num_particles, num_dimensions))
velocities = np.zeros((num_particles, num_dimensions))
p_best = particles.copy()
p_best_scores = np.full(num_particles, np.inf)
g_best = particles[0].copy()
g_best_score = np.inf

# PSO 主循环
for _ in range(max_iter):
    for i in range(num_particles):
        current_score = float(np.asarray(fitness(particles[i, :])).item())
        if current_score < p_best_scores[i]:
            p_best_scores[i] = current_score
            p_best[i, :] = particles[i, :]
        if current_score < g_best_score:
            g_best_score = current_score
            g_best = particles[i, :].copy()

    for i in range(num_particles):
        velocities[i, :] = (
            w * velocities[i, :]
            + c1 * rng.random() * (p_best[i, :] - particles[i, :])
            + c2 * rng.random() * (g_best - particles[i, :])
        )
        velocities[i, :] = np.maximum(np.minimum(velocities[i, :], ub), lb)
        particles[i, :] = particles[i, :] + velocities[i, :]
        particles[i, :] = np.maximum(np.minimum(particles[i, :], ub), lb)

    if float(np.std(p_best_scores)) < tolerance:
        break

V_mppt = float(g_best[0])
I_mppt = float(pv_iv(V_mppt))
P_mppt = V_mppt * I_mppt
print(f"PSO找到的最优电压: {V_mppt:.2f} V")
print(f"对应电流: {I_mppt:.2f} A, 最大功率: {P_mppt:.2f} W")