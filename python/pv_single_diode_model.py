"""
pv_single_diode_model.py — 光伏电池单二极管模型，绘制 I-V 与 P-V 曲线
对应 pv_single_diode_model.m：计算并绘制 I-V/P-V 特性，输出最大功率点。
运行：uv run python pv_single_diode_model.py
"""
import numpy as np
import plot

from pv_model import pv_params, pv_current

pv = pv_params()
V_oc = pv["Voc_stc"]

V_PV = np.linspace(0, V_oc, 1000)
I_PV = pv_current(V_PV, pv)
P = V_PV * I_PV

i_max = int(np.argmax(P))
V_mppt, P_mppt = V_PV[i_max], P[i_max]

print(f"开路电压: {V_oc:.2f} V, 短路电流: {pv['Isc_stc']:.2f} A")
print(f"最大功率点: V={V_mppt:.2f} V, P={P_mppt:.2f} W")

# I-V 特性曲线
fig, ax = plot.new_figure("pv_single_diode_iv")
ax.plot(V_PV, I_PV, "-b", label="I-V")
ax.set_xlabel("Voltage (V)")
ax.set_ylabel("Current (A)")
ax.set_title("I-V 特性曲线")
ax.legend()
plot.save(fig)

# P-V 特性曲线
fig, ax = plot.new_figure("pv_single_diode_pv")
ax.plot(V_PV, P, "-r", label="P-V")
ax.plot(V_mppt, P_mppt, "ko", markersize=8, markerfacecolor="k", label="MPP")
ax.set_xlabel("Voltage (V)")
ax.set_ylabel("Power (W)")
ax.set_title("P-V 特性曲线")
ax.legend()
plot.save(fig)