"""
pv_iv_family.py — 光伏 I-V/P-V 曲线族（不同辐照度、不同温度）
对应理论文档 1.2 节"环境因素对特性曲线的影响"的配图。
模型统一来自 pv_model.pv_params() + pv_current()。
运行：uv run python pv_iv_family.py
"""
import numpy as np
import plot

from pv_model import pv_params, pv_current

pv = pv_params()
V = np.linspace(0, pv["Voc_stc"], 500)

# ---- 不同辐照度下的 I-V / P-V 曲线族 ----
fig, (ax1, ax2) = plot.new_subplots("pv_irradiance_family", 1, 2, figsize=(11, 4))
for G in [200, 400, 600, 800, 1000]:
    I = pv_current(V, pv, G, 25)
    P = V * I
    ax1.plot(V, I, label=f"G={G} W/m²")
    ax2.plot(V, P, label=f"G={G} W/m²")
ax1.set_xlabel("Voltage (V)")
ax1.set_ylabel("Current (A)")
ax1.set_title("不同辐照度下 I-V")
ax1.grid(True)
ax1.legend()
ax2.set_xlabel("Voltage (V)")
ax2.set_ylabel("Power (W)")
ax2.set_title("不同辐照度下 P-V")
ax2.grid(True)
ax2.legend()
plot.save(fig)

# ---- 不同温度下的 I-V / P-V 曲线族 ----
fig, (ax1, ax2) = plot.new_subplots("pv_temperature_family", 1, 2, figsize=(11, 4))
for T in [0, 25, 50, 75]:
    I = pv_current(V, pv, 1000, T)
    P = V * I
    ax1.plot(V, I, label=f"T={T}°C")
    ax2.plot(V, P, label=f"T={T}°C")
ax1.set_xlabel("Voltage (V)")
ax1.set_ylabel("Current (A)")
ax1.set_title("不同温度下 I-V")
ax1.grid(True)
ax1.legend()
ax2.set_xlabel("Voltage (V)")
ax2.set_ylabel("Power (W)")
ax2.set_title("不同温度下 P-V")
ax2.grid(True)
ax2.legend()
plot.save(fig)

print("已生成辐照度族与温度族曲线图（figures/pv_irradiance_family.png, pv_temperature_family.png）")