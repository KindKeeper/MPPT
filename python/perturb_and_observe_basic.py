"""
perturb_and_observe_basic.py — 扰动观察法（P&O）— 基于特性曲线表格实现
对应 perturb_and_observe_basic.m：先生成单二极管模型 I-V/P-V 曲线，
再在曲线上用变步长 P&O 按索引搜索最大功率点。
运行：uv run python perturb_and_observe_basic.py
"""
import numpy as np
import plot

from pv_model import pv_params, pv_current

pv = pv_params()
V_oc = pv["Voc_stc"]
P_threshold = 2   # 功率阈值(W)
V_threshold = 1   # 电压阈值(V)

# 变步长系数
k1 = 0.2          # 远离峰值时大步长
k2 = 0.02         # 接近峰值时小步长

# 光伏特性曲线（统一单二极管模型）
V_PV = np.linspace(0, V_oc, 1000)
I_PV = np.maximum(pv_current(V_PV, pv), 0)
P = V_PV * I_PV
V_vec, I_vec, P_vec = V_PV, I_PV, P

# 算法初始化
N = V_vec.size
dV = float(np.mean(np.diff(V_vec)))
iter_ = 1000
idx = 0            # Python 从 0 开始索引（对应 .m 中 1 起步，整体左移 1）
V_s = np.zeros(iter_)
P_s = np.zeros(iter_)

# 变步长 P&O 主循环
for k in range(iter_):
    V = V_vec[idx]
    P = P_vec[idx]
    V_new = V_vec[idx + 1]
    P_new = P_vec[idx + 1]
    dP = P_new - P
    dV = V_new - V

    V_s[k] = V
    P_s[k] = P

    # 变步长策略
    if abs(dP) > P_threshold:
        step = k1 * abs(dP)
    elif abs(dV) < V_threshold:
        step = k2 * abs(dP)
    else:
        step = 0.1

    step_idx = max(1, round(abs(step) / max(abs(dV), 1e-9)))
    step_idx = int(np.sign(P_new - P)) * step_idx
    idx_new = idx + step_idx
    idx = min(max(idx_new, 0), N - 2)   # 上限 N-2 对应 .m 的 N-1

V_mppt = V_vec[idx]
P_mppt = P_vec[idx]

print(f"MPPT电压: {V_mppt:.2f} V")
print(f"MPPT功率: {P_mppt:.2f} W")

# 可视化
fig, ax = plot.new_figure("perturb_po_iv")
ax.plot(V_vec, I_vec, "-b", label="I-V")
ax.plot(V_mppt, I_vec[idx], "ro", label="MPP")
ax.set_xlabel("Voltage (V)")
ax.set_ylabel("Current (A)")
ax.set_title("I–V")
ax.legend()
plot.save(fig)

fig, ax = plot.new_figure("perturb_po_pv")
ax.plot(V_vec, P_vec, "-r", label="P-V")
ax.plot(V_mppt, P_mppt, "ko", label="MPP")
ax.set_xlabel("Voltage (V)")
ax.set_ylabel("Power (W)")
ax.set_title("P–V")
ax.legend()
plot.save(fig)

fig, ax = plot.new_figure("perturb_po_p_iter")
ax.plot(P_s, linewidth=1.2)
ax.set_xlabel("迭代过程")
ax.set_ylabel("P (W)")
plot.save(fig)

fig, ax = plot.new_figure("perturb_po_v_iter")
ax.plot(V_s, linewidth=1.2)
ax.set_xlabel("迭代过程")
ax.set_ylabel("$V_{pv}$ (V)")
plot.save(fig)