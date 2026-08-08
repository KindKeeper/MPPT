"""
plot.py — 统一的出图助手（将曲线保存到 python/figures/ 而非弹窗）
用于无界面环境（服务器/CI），保证所有脚本图和 .m 版本一一对应。
"""
import os
import matplotlib

matplotlib.use("Agg")  # 无界面后端

import matplotlib.pyplot as plt

# 配置中文字体（Noto Sans CJK），避免图形中文显示为方块
import matplotlib.font_manager as fm

_CJK_FONTS = [
    "Noto Sans CJK SC",
    "Noto Sans CJK TC",
    "Noto Serif CJK SC",
    "WenQuanYi Micro Hei",
    "Source Han Sans SC",
    "AR PL UMing CN",
]
_candidates = {f.name for f in fm.fontManager.ttflist}
for _f in _CJK_FONTS:
    if _f in _candidates:
        plt.rcParams["font.sans-serif"] = [_f] + list(plt.rcParams["font.sans-serif"])
        plt.rcParams["axes.unicode_minus"] = False  # 正常显示负号
        break

_FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")


def new_figure(name):
    """创建图形窗口，返回 (fig, ax)。"""
    os.makedirs(_FIG_DIR, exist_ok=True)
    fig, ax = plt.subplots()
    ax.grid(True)
    fig._mpptsave = _save_path(name)
    return fig, ax


def new_subplots(name, nrows, ncols, figsize=None):
    """创建多子图图形窗口，返回 (fig, axes)。"""
    os.makedirs(_FIG_DIR, exist_ok=True)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    fig._mpptsave = _save_path(name)
    return fig, axes


def _save_path(name):
    return os.path.join(_FIG_DIR, name if name.endswith(".png") else name + ".png")


def save(fig):
    """保存并关闭图形到 figures/ 目录。"""
    os.makedirs(_FIG_DIR, exist_ok=True)
    path = fig._mpptsave
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"已保存图形: {path}")
    return path