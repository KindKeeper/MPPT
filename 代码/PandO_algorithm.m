%% 扰动观察法（P&O）— 基础单二极管模型版
% 功能：通过固定/变步长扰动电压，观察功率变化方向，逐步逼近最大功率点
% 运行环境：MATLAB 或 Octave
% 说明：本脚本为独立完整示例，直接运行即可看到搜索结果。
%       光伏模型参数统一来自 pv_params() + pv_current()。

clear; clc; close all;

% ===== 光伏组件参数（统一参数源） =====
pv = pv_params();
V_oc = pv.Voc_stc;

% 单二极管模型电流-电压关系（统一模型入口）
pv_iv = @(V) pv_current(V, pv);

% ===== P&O 参数 =====
step = 0.05;        % 初始扰动步长 [V]
P_threshold = 1;    % 功率变化阈值 [W]，超过则用大步长
tol = 1e-3;         % 收敛容差 [W]
max_iter = 1000;    % 最大迭代次数

% 初始化：从约 1/2 开路电压处开始搜索
V = V_oc/2;
I = pv_iv(V);
P = V * I;

% P&O 主循环：固定步长扰动，依据功率增减反转方向（含限幅）
for iter = 1:max_iter
    V_new = V + step;                 % 施加扰动
    I_new = pv_iv(V_new);
    P_new = V_new * I_new;

    if P_new > P
        step = abs(step);             % 功率上升，保持方向
    else
        step = -abs(step);            % 功率下降，反转方向
    end
    V = min(max(V_new, 0), V_oc);        % 限幅

    if abs(P_new - P) < tol
        break;                           % 收敛
    end
    P = P_new;
end

V_mppt = V; I_mppt = pv_iv(V); P_mppt = V * I_mppt;

fprintf('MPPT电压: %.2f V\n', V_mppt);
fprintf('MPPT电流: %.2f A\n', I_mppt);
fprintf('MPPT功率: %.2f W\n', P_mppt);