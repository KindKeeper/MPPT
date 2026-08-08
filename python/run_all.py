#!/usr/bin/env python
"""
run_all.py — MPPT 仓库 Python 版一键运行脚本
功能：依次运行 python/ 目录下所有可独立运行的 MPPT 算法脚本，验证可复现。
用法：在仓库根目录执行  uv run python python/run_all.py
说明：每个脚本以独立子进程运行，互不干扰。
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))

runnable = [
    "pv_single_diode_model",      # 1. 光伏单二极管模型与 I-V/P-V 曲线
    "pv_iv_family",               # 2. 不同辐照度/温度下的曲线族
    "inc_algorithm",              # 3. 电导增量法
    "pando_algorithm",            # 4. 扰动观察法（基础版）
    "perturb_and_observe_basic",  # 5. 扰动观察法（曲线实现）
    "global_scan_pv",             # 6. 全局扫描+精跟（局部遮阴多峰）
    "mppt_adaptive_slope_boost",  # 7. 变步长自适应 P&O + Boost（核心）
    "pso_algorithm",              # 8. 粒子群算法
    "mppt_algorithm_comparison",  # 9. P&O vs PSO 对比
]

print("========== MPPT 仓库(Python)一键验证开始 ==========")
n_pass = n_fail = 0

for i, script in enumerate(runnable, 1):
    fpath = os.path.join(HERE, script + ".py")
    if not os.path.isfile(fpath):
        print(f"跳过：找不到文件 {fpath}")
        continue
    print(f"\n----- [{i}/{len(runnable)}] 运行 {script} -----")
    logf = os.path.join(tempfile.gettempdir(), f"run_{script}.log")
    with open(logf, "w") as f:
        proc = subprocess.run(
            [sys.executable, fpath],
            cwd=HERE,
            stdout=f,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if proc.returncode == 0:
        n_pass += 1
        print(f"     ✓ {script} 运行完成")
    else:
        n_fail += 1
        print(f"     ✗ {script} 运行失败（exit={proc.returncode}）")
        with open(logf) as f:
            for line in f.readlines()[:6]:
                print(f"        | {line.rstrip()}")

print(f"\n========== 验证结束：通过 {n_pass} / {n_pass + n_fail}，失败 {n_fail} ==========")
