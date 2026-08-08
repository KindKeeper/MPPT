"""
mppt_algorithm_comparison.py — 光伏 MPPT 算法性能对比：P&O vs PSO
对应 mppt_algorithm_comparison.m：在统一单二极管模型下分别运行
扰动观察法(P&O)和粒子群(PSO)寻找最大功率点，对比收敛结果与迭代开销。
运行：uv run python mppt_algorithm_comparison.py
"""
import numpy as np
import plot

from pv_model import pv_params, pv_current

rng = np.random.default_rng(42)

pv = pv_params()
V_oc = pv["Voc_stc"]
pv_eq = lambda V: pv_current(V, pv)

# ---- 方法一：扰动观察法(P&O) ----
delta = 0.05
V = 1.0
P_po_hist = []
for _ in range(5000):
    I = pv_eq(V)
    P = V * I
    P_po_hist.append(P)
    V_new = V + delta
    I_new = pv_eq(V_new)
    P_new = V_new * I_new
    if P_new > P:
        V = V_new
    else:
        V = V - delta
    V = min(max(V, 0), V_oc)
    if abs(P_new - P) < 1e-3:
        break
V_po = V
I_po = pv_eq(V_po)
P_po = V_po * I_po

# ---- 方法二：粒子群算法(PSO) ----
num_particles, max_iter = 30, 50
w, c1, c2 = 0.5, 1.5, 2.0
lb, ub = 0, 36.6
particles = lb + (ub - lb) * rng.random(num_particles)
velocities = np.zeros(num_particles)
p_best = particles.copy()
p_best_score = np.array([v * pv_eq(v) for v in particles])
gi = int(np.argmax(p_best_score))
g_best = particles[gi]
g_best_score = p_best_score[gi]
for _ in range(max_iter):
    for i in range(num_particles):
        score = particles[i] * pv_eq(particles[i])
        if score > p_best_score[i]:
            p_best_score[i] = score
            p_best[i] = particles[i]
        if score > g_best_score:
            g_best_score = score
            g_best = particles[i]
    for i in range(num_particles):
        r1, r2 = rng.random(), rng.random()
        velocities[i] = (
            w * velocities[i]
            + c1 * r1 * (p_best[i] - particles[i])
            + c2 * r2 * (g_best - particles[i])
        )
        particles[i] = max(lb, min(ub, particles[i] + velocities[i]))
V_pso = g_best
I_pso = pv_eq(V_pso)
P_pso = V_pso * I_pso

# ---- 结果展示与对比 ----
V_sweep = np.linspace(0, V_oc, 500)
P_sweep = V_sweep * pv_eq(V_sweep)

fig, ax = plot.new_figure("comparison")
ax.plot(V_sweep, P_sweep, "-", linewidth=1.5, label="P-V曲线")
ax.plot(V_po, P_po, "bo", markersize=8, markerfacecolor="b", label="P&O 结果")
ax.plot(V_pso, P_pso, "rs", markersize=8, markerfacecolor="r", label="PSO 结果")
ax.set_xlabel("V (V)")
ax.set_ylabel("P (W)")
ax.set_title("MPPT 算法对比")
ax.legend()
plot.save(fig)

print("—— 对比结果 ——")
print(f"P&O : V={V_po:.2f}V, I={I_po:.2f}A, P={P_po:.2f}W")
print(f"PSO : V={V_pso:.2f}V, I={I_pso:.2f}A, P={P_pso:.2f}W")
print(f"P&O 迭代次数: {len(P_po_hist)}")
print(f"PSO 迭代代数: {max_iter}")