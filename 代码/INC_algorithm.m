%% 电导增量法（INC）— 简化教学版
% 功能：通过电导增量法寻找光伏电池的最大功率点（MPPT）
% 原理：dP/dV = 0 处为最大功率点，用增量电导近似判定
% 运行环境：MATLAB 或 Octave
% 说明：本脚本刻意使用线性简化模型（I=Isc-(Isc/Voc)*V），
%       以便清晰展示 INC 算法"dP/dV=0 收敛"的核心思想，
%       经典结论：线性模型下最大功率点恰在 Voc/2 处。
%       使用完整单二极管模型（考虑辐照/温度）的实现见
%       mppt_algorithm_comparison.m（统一 pv_params/pv_current 模型）。

clear; clc; close all;

% ===== 光伏电池参数（简化线性模型） =====
Voc = 36;        % 开路电压 [V]
Isc = 5;         % 短路电流 [A]
delta = 0.01;    % 扰动步长 [V]
tolerance = 1e-6;% 收敛容差
max_iter = 1000; % 最大迭代次数（防止极端情况死循环）

% 初始电压猜测（取开路电压一半）
Vpv = Voc/2;

% INC 主循环
for iter = 1:max_iter
    % 简化电流模型：电流随电压线性下降
    Ipv = Isc - (Isc/Voc)*Vpv;
    % 电流对电压的导数
    dI_dV = -Isc/Voc;
    % 功率对电压的导数：dP/dV = I + V*dI/dV
    dP_dV = Ipv + Vpv*dI_dV;

    % 依据 dP/dV 符号调整电压方向
    if dP_dV > tolerance
        Vpv = Vpv + delta;   % 沿功率上升方向调整
    elseif dP_dV < -tolerance
        Vpv = Vpv - delta;   % 反向调整
    else
        break;               % 收敛到 MPP
    end

    % 电压限幅在有效范围 [0, Voc]
    Vpv = min(max(Vpv, 0), Voc);
end

fprintf('最大功率点电压: %.2f V\n', Vpv);
fprintf('对应电流: %.2f A, 功率: %.2f W\n', ...
        Isc - (Isc/Voc)*Vpv, Vpv*(Isc - (Isc/Voc)*Vpv));