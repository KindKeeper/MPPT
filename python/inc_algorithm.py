"""
inc_algorithm.py — 电导增量法（INC）— 简化教学版
对应 INC_algorithm.m：通过 dP/dV=0 收敛寻找最大功率点。
刻意使用线性简化模型（I=Isc-(Isc/Voc)*V），经典结论：MPP 恰在 Voc/2 处。
运行：uv run python inc_algorithm.py
"""
Voc = 36        # 开路电压 [V]
Isc = 5         # 短路电流 [A]
delta = 0.01    # 扰动步长 [V]
tolerance = 1e-6
max_iter = 1000

Vpv = Voc / 2   # 初始电压猜测（开路电压一半）

# INC 主循环
for _ in range(max_iter):
    Ipv = Isc - (Isc / Voc) * Vpv      # 简化电流模型
    dI_dV = -Isc / Voc                 # 电流对电压的导数
    dP_dV = Ipv + Vpv * dI_dV         # dP/dV = I + V*dI/dV

    if dP_dV > tolerance:
        Vpv += delta
    elif dP_dV < -tolerance:
        Vpv -= delta
    else:
        break                        # 收敛到 MPP

    Vpv = min(max(Vpv, 0), Voc)      # 电压限幅

Ipv_f = Isc - (Isc / Voc) * Vpv
print(f"最大功率点电压: {Vpv:.2f} V")
print(f"对应电流: {Ipv_f:.2f} A, 功率: {Vpv * Ipv_f:.2f} W")