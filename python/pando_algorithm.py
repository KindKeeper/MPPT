"""
pando_algorithm.py — 扰动观察法（P&O）— 基础单二极管模型版
对应 PandO_algorithm.m：通过扰动电压观察功率变化方向，逐步逼近最大功率点。
模型统一来自 pv_model.pv_params() + pv_current()。
运行：uv run python pando_algorithm.py
"""
from pv_model import pv_params, pv_current

pv = pv_params()
V_oc = pv["Voc_stc"]
pv_iv = lambda V: pv_current(V, pv)

# P&O 参数
step = 0.05        # 初始扰动步长 [V]
tol = 1e-3         # 收敛容差 [W]
max_iter = 1000

# 初始化：从约 1/2 开路电压处开始
V = V_oc / 2
I = pv_iv(V)
P = V * I

# P&O 主循环
for _ in range(max_iter):
    V_new = V + step
    I_new = pv_iv(V_new)
    P_new = V_new * I_new

    if P_new > P:
        step = abs(step)       # 功率上升，保持方向
    else:
        step = -abs(step)      # 功率下降，反转方向
    V = min(max(V_new, 0), V_oc)

    if abs(P_new - P) < tol:
        break
    P = P_new

I_mppt = pv_iv(V)
P_mppt = V * I_mppt
print(f"MPPT电压: {V:.2f} V")
print(f"MPPT电流: {I_mppt:.2f} A")
print(f"MPPT功率: {P_mppt:.2f} W")