"""
pv_model.py — 光伏单二极管模型统一核心（对应 MATLAB/Octave 的 pv_params + pv_current + 辅助函数）
功能：全仓库唯一的 PV 物理模型来源。所有算法脚本只依赖本模块，保证跨语言结果可比。
运行环境：Python 3 (numpy)

对应 .m 文件：
    pv_params.m                -> pv_params()
    pv_current.m               -> pv_current()
    pv_current_given_voltage.m -> pv_current_given_voltage()
    solve_pv_operating_point.m -> solve_pv_operating_point()
    sweep_pv_iv.m              -> sweep_pv_iv()
"""
import numpy as np


def pv_params(name=None):
    """返回标准 PV 组件物理参数（与 pv_params.m 完全一致）。

    参数
    ----
    name : str, optional
        预留：未来支持多型号组件。当前仅默认 60 片单晶硅组件。

    返回
    ----
    pv : dict
        含 Ns, Isc_stc, Voc_stc, n, Rs, Rsh, alpha_Isc, Eg, q, k 等字段。
    """
    pv = {
        "Ns": 60,          # 串联电池片数
        "Isc_stc": 8.7,    # STC 短路电流 [A]
        "Voc_stc": 37.0,   # STC 开路电压 [V]
        "n": 1.3,          # 二极管理想因子
        "Rs": 0.35,        # 串联电阻 [Ohm]
        "Rsh": 500,        # 并联电阻 [Ohm]
        "alpha_Isc": 0.005,  # 电流温度系数 [A/°C]
        "Eg": 1.12,        # 硅禁带宽度 [eV]
        "q": 1.602176634e-19,  # 元电荷 [C]
        "k": 1.380649e-23,     # 玻尔兹曼常数 [J/K]
    }
    return pv


def pv_current(V, pv, G=1000, T=25):
    """给定电压，用单二极管模型计算光伏组件输出电流（统一模型入口）。

    参数
    ----
    V : float or np.ndarray
        电压 [V]（标量或向量）
    pv : dict
        由 pv_params() 生成的参数
    G : float
        辐照度 [W/m^2]，默认 1000（STC）
    T : float
        组件温度 [°C]，默认 25（STC）

    返回
    ----
    I : float or np.ndarray
        电流 [A]，与 V 同形状，物理上截断为非负
    """
    V = np.maximum(V, 0)

    # 热电压（考虑串联片数）
    Vt = pv["Ns"] * pv["k"] * (T + 273.15) / pv["q"]

    # 光生电流：随辐照度线性变化 + 温度微调
    Iph = pv["Isc_stc"] * (G / 1000) + pv["alpha_Isc"] * (T - 25)

    # 反向饱和电流（由 STC 条件反推）
    Vt_stc = pv["Ns"] * pv["k"] * (25 + 273.15) / pv["q"]
    I0_stc = pv["Isc_stc"] / (np.exp(pv["Voc_stc"] / (pv["n"] * Vt_stc)) - 1)

    # I0 温度修正（含禁带宽度 Eg 影响）
    T_K = T + 273.15
    I0 = I0_stc * (T_K / 298.15) ** 3 * np.exp(
        (pv["q"] * pv["Eg"] / pv["k"]) * (1 / 298.15 - 1 / T_K) / (pv["n"] * pv["Ns"])
    )

    # 单二极管隐式方程（固定点迭代求解，对典型组件收敛快速、数值稳定）
    # I = Iph - I0*(exp((V + I*Rs)/(n*Vt)) - 1) - (V + I*Rs)/Rsh
    I = Iph
    for _ in range(30):
        I = Iph - I0 * (np.exp((V + I * pv["Rs"]) / (pv["n"] * Vt)) - 1) - (
            V + I * pv["Rs"]
        ) / pv["Rsh"]
        I = np.maximum(I, 0)
    return I


def pv_current_given_voltage(V, G, T, params):
    """单二极管模型电流计算的薄封装（兼容 mppt_adaptive_slope_boost 接口）。

    对应 pv_current_given_voltage.m，内部调用 pv_current。
    """
    return pv_current(V, params, G, T)


def _bisection(f, a, b, tol=1e-10, max_iter=200):
    """在 [a, b] 上二分法求 f(x)=0 的根（要求 f(a)、f(b) 异号）。

    对应 solve_pv_operating_point 中 fzero 的区间搜索模式，
    手写二分避免引入 scipy 依赖，保持项目轻量。
    """
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        raise ValueError("二分法要求区间端点函数值异号")
    for _ in range(max_iter):
        c = (a + b) / 2
        fc = f(c)
        if abs(fc) < tol or (b - a) / 2 < tol:
            return c
        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return (a + b) / 2


def solve_pv_operating_point(D, Rload, G, T, params):
    """求解给定占空比下的光伏工作点电压和电流。

    对应 solve_pv_operating_point.m：
    理想 Boost 变换器在光伏侧看到的等效输入电阻 Rin = Rload * (1-D)^2，
    工作点满足 f(V) = Ipv(V) - V/Rin = 0。

    返回 (Vsol, Isol)。
    """
    Rin = Rload * (1 - D) ** 2
    if Rin <= 0:
        Rin = 1e-6  # 防零除保护

    f = lambda V: pv_current_given_voltage(max(V, 0), G, T, params) - max(V, 0) / Rin

    # 简化开路电压估计（与 .m 一致，暂不做温度修正）
    Voc_est = params["Voc_stc"]
    V0 = 0.7 * Voc_est  # 初始猜测

    # 区间搜索：f(0)=Ipv(0)>0，f(1.2*Voc_est)≈0-正数<0，端点异号
    a, b = 0.0, max(1e-6, 1.2 * Voc_est)
    try:
        Vsol = _bisection(f, a, b)
    except ValueError:
        # 区间搜索失败时回退：围绕单点猜测微扰二分
        a, b = max(0.0, V0 * 0.5), min(1.2 * Voc_est, V0 * 1.5)
        if f(a) * f(b) > 0:
            Vsol = V0  # 仍失败则返回猜测值
        else:
            Vsol = _bisection(f, a, b)
    Vsol = max(Vsol, 0)  # 确保电压非负
    Isol = pv_current_given_voltage(Vsol, G, T, params)
    return Vsol, Isol


def sweep_pv_iv(G, T, params):
    """扫描生成完整 I-V 特性曲线。

    对应 sweep_pv_iv.m：电压从 0 到 1.15*Voc_stc 均匀取 200 点。
    返回 (Vvec, Ivec)。
    """
    Vmax = params["Voc_stc"] * 1.15
    Vvec = np.linspace(0, Vmax, 200)
    Ivec = np.array([pv_current_given_voltage(v, G, T, params) for v in Vvec])
    return Vvec, Ivec
