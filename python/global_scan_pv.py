"""
global_scan_pv.py — 全局扫描 + 精跟 MPPT（应对局部遮阴多峰）
对应 global_scan_pv.m：模拟局部遮阴下 P-V 多峰，对比"全局扫描+精调"
与"普通 P&O（陷入局部峰）"两种策略。
用"两串并联"教学模型：各串辐照度不同、各有自身峰值，叠加形成双峰。
运行：uv run python global_scan_pv.py
"""
import numpy as np
import plot

from pv_model import pv_params, pv_current

pv = pv_params()
V_oc = pv["Voc_stc"]
pv_iv = lambda V, G: pv_current(V, pv, G, 25)

# 构建"局部遮阴"双峰 P-V 曲线
G_active = 1000    # 光照充足一串
G_shaded = 300     # 遮阴一串
N = 1000
V_sweep = np.linspace(0, V_oc, N)

I_hi = pv_iv(V_sweep, G_active)
I_lo = pv_iv(V_sweep, G_shaded)
I_partial = 1.0 * I_hi + 0.7 * I_lo        # 光照串为主，遮阴串配角
P_partial = V_sweep * I_partial            # 双峰功率曲线

# 线性插值辅助（对应 MATLAB interp1(...,'linear',0)）
def interp_lin(x_new):
    return np.interp(x_new, V_sweep, P_partial, left=0.0, right=0.0)

# ---- 纯局部 P&O：从高压侧出发，可能陷入局部峰值 ----
V0 = 36.0
delta = -0.5
V_po = V0
for _ in range(500):
    P0 = interp_lin(V_po)
    V1 = V_po + delta
    P1 = interp_lin(V1)
    if P1 > P0:
        V_po = V1
    else:
        V_po = V1 - delta
        delta = -delta
    V_po = min(max(V_po, 0), V_oc)
P_po_trap = interp_lin(V_po)

# ---- 全局扫描策略：先扫全范围锁定全局峰，再精调 ----
i_global = int(np.argmax(P_partial))
V_global, P_global = V_sweep[i_global], P_partial[i_global]

V_fine = np.linspace(max(0, V_global - 0.5), min(V_oc, V_global + 0.5), 200)
P_fine = np.interp(V_fine, V_sweep, P_partial, left=0.0, right=0.0)
i_fine = int(np.argmax(P_fine))
V_mppt, P_mppt = V_fine[i_fine], P_fine[i_fine]

# ---- 结果对比 ----
print("===== 局部遮阴下 MPPT 结果对比 =====")
print(f"局部P&O(易陷局部峰): V={V_po:.2f} V, P={P_po_trap:.2f} W")
print(f"全局扫描找全局峰:   V={V_global:.2f} V, P={P_global:.2f} W")
print(f"扫描+精调最终MPP:    V={V_mppt:.2f} V, P={P_mppt:.2f} W")
print("结论：全局扫描策略优于局部P&O，避免陷入局部峰值。")

# ---- 可视化 ----
fig, ax = plot.new_figure("global_scan_pv")
ax.plot(V_sweep, P_partial, "-b", linewidth=1.5, label="多峰P-V曲线")
ax.plot(V_mppt, P_mppt, "ro", markersize=8, markerfacecolor="r", label="全局扫描MPPT")
ax.plot(V_po, P_po_trap, "ks", markersize=8, markerfacecolor="k", label="局部P&O(陷阱)")
ax.set_xlabel("Voltage (V)")
ax.set_ylabel("Power (W)")
ax.set_title("局部遮阴下 MPPT：全局扫描 vs 局部P&O")
ax.legend()
plot.save(fig)