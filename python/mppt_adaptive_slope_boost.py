"""
mppt_adaptive_slope_boost.py — PV 单二极管模型 + 理想 Boost 平均模型
                              + 变步长（基于 dP/dV）P&O MPPT
对应 mppt_adaptive_slope_boost.m，依赖 pv_model 提供的辅助函数。
运行：uv run python mppt_adaptive_slope_boost.py
"""
import numpy as np
import plot

from pv_model import (
    pv_params,
    solve_pv_operating_point,
    sweep_pv_iv,
)

pv = pv_params()

# 环境参数
G = 800       # 辐照度 [W/m^2]
Tc = 35       # 组件温度 [°C]

# 负载与 Boost 平均模型
Rload = 20            # 直流负载电阻 [Ohm]
Dmin, Dmax = 0.02, 0.95
D = 0.5

# ---- MPPT（基于斜率 dP/dV 的自适应扰动）参数 ----
Nstep = 300
Kgain = 1e-3
step_min, step_max = 2e-4, 1.5e-2
slope_eps = 1e-3
slope_LP = 0.7

Dlog = np.zeros(Nstep)
Vlog = np.zeros(Nstep)
Ilog = np.zeros(Nstep)
Plog = np.zeros(Nstep)

Vk, Ik = solve_pv_operating_point(D, Rload, G, Tc, pv)
Pk = Vk * Ik
dPdV_est_f = 0.0

# ---- MPPT 主循环 ----
for k in range(Nstep):
    Dlog[k], Vlog[k], Ilog[k], Plog[k] = D, Vk, Ik, Pk

    D_probe = min(max(D + 1e-3, Dmin), Dmax)
    V_probe, I_probe = solve_pv_operating_point(D_probe, Rload, G, Tc, pv)
    P_probe = V_probe * I_probe

    dV = V_probe - Vk
    dP = P_probe - Pk
    dPdV = (abs(dV) > 1e-9) * (dP / (dV + (abs(dV) <= 1e-9)))
    dPdV_est_f = slope_LP * dPdV_est_f + (1 - slope_LP) * dPdV

    if abs(dPdV_est_f) < slope_eps:
        dir_ = 0.0                   # 接近 MPP，方向为零
        st = step_min                # 最小步长，降低振荡
    else:
        dir_ = np.sign(dPdV_est_f)   # 斜率决定电压调整方向
        st = min(step_max, max(step_min, Kgain * abs(dPdV_est_f)))

    D = D - dir_ * st
    D = min(max(D, Dmin), Dmax)

    Vk, Ik = solve_pv_operating_point(D, Rload, G, Tc, pv)
    Pk = Vk * Ik

# ---- 结果可视化 ----
Vvec, Ivec = sweep_pv_iv(G, Tc, pv)
Pvec = Vvec * Ivec

fig, ax = plot.new_figure("boost_iv")
ax.plot(Vvec, Ivec, linewidth=1.5, label="I–V")
ax.plot(Vlog[-1], Ilog[-1], "o", label="MPPT point")
ax.set_xlabel("$V_{pv}$ (V)")
ax.set_ylabel("$I_{pv}$ (A)")
ax.set_title("PV I–V (current G,T) & MPPT result")
ax.legend()
plot.save(fig)

fig, ax = plot.new_figure("boost_pv")
ax.plot(Vvec, Pvec, linewidth=1.5, label="P–V")
ax.plot(Vlog[-1], Plog[-1], "o", label="MPPT point")
ax.set_xlabel("$V_{pv}$ (V)")
ax.set_ylabel("$P_{pv}$ (W)")
ax.set_title("PV P–V (current G,T) & MPPT result")
ax.legend()
plot.save(fig)

fig, ax = plot.new_figure("boost_power")
ax.plot(Plog, linewidth=1.2)
ax.set_xlabel("Iteration")
ax.set_ylabel("P (W)")
ax.set_title("Convergence of Power")
plot.save(fig)

fig, ax = plot.new_figure("boost_voltage")
ax.plot(Vlog, linewidth=1.2)
ax.set_xlabel("Iteration")
ax.set_ylabel("$V_{pv}$ (V)")
ax.set_title("PV Voltage vs. Iteration")
plot.save(fig)

fig, ax = plot.new_figure("boost_duty")
ax.plot(Dlog, linewidth=1.2)
ax.set_xlabel("Iteration")
ax.set_ylabel("D")
ax.set_title("Duty Ratio vs. Iteration")
plot.save(fig)

print(f"MPPT结果: V={Vlog[-1]:.2f} V, P={Plog[-1]:.2f} W")