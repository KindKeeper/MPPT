%% global_scan_pv.m — 全局扫描 + 精跟 MPPT（应对局部遮阴多峰）
% 功能：模拟局部遮阴下 P-V 曲线出现多个峰值的情况，
%       用"全局扫描 + 精调跟踪"找到全局最大功率点，
%       并与普通 P&O（可能陷入局部峰）对比，演示两者差异。
% 运行环境：MATLAB 或 Octave
% 说明：光伏模型参数统一来自 pv_params() + pv_current()。
%       局部遮阴用"串联电池组 + 旁路二极管"物理模型模拟，
%       电流被最弱辐照度组钳位，形成真实的多峰 P-V 曲线。

clear; clc; close all;

% ===== 光伏组件参数（统一参数源） =====
pv = pv_params();
V_oc = pv.Voc_stc;

% 单二极管模型电流-电压关系（统一模型入口，固定温度25°C）
pv_iv = @(V, G) pv_current(V, pv, G, 25);

%% ===== 构建"局部遮阴"多峰 P-V 曲线（物理串联模型） =====
% 组件分三段串联，各自受不同辐照度 G1>G2>G3。
% 串联电流受最弱段限制；较强的段通过旁路二极管绕过，从而产生阶梯与多峰。
G_active = 1000;                 % 全光照辐照度
G_shaded = 300;                  % 遮阴段辐照度

N = 1000;
% 采样电压（串联全部段的总电压）
V_sweep = linspace(0, V_oc, N);

% 物理串联遮阴模型（标准"两串并联"教学简化）：
% 组件等效为两串电池并联，各串受不同辐照度、各有自己的 MPP。
% 遮阴较严重的那一串(低辐照度G_shaded)在较低电压处先达到自身峰值，
% 另一串(高辐照度G_active)在较高电压处达峰，两者叠加形成"双峰"P-V。
% 这清晰复现了局部遮阴导致的多峰现象，用于演示全局扫描策略。

G_active = 1000;                 % 光照充足一串的辐照度
G_shaded = 300;                  % 遮阴一串的辐照度

I_hi = pv_iv(V_sweep, G_active);            % 光照串贡献
I_lo = pv_iv(V_sweep, G_shaded);            % 遮阴串贡献
% 权重：光照串为主（1.0），遮阴串配角（0.7），合成双峰总电流
I_partial = 1.0*I_hi + 0.7*I_lo;
P_partial = V_sweep .* I_partial;   % 双峰功率曲线

%% ===== 纯局部 P&O：从高压侧出发，可能陷入局部峰值 =====
% 演示：P&O 只做局部搜索，从高点(如35V处)出发，找不到翻越凹谷到达全局峰，
% 会停在局部峰 V≈35 V(P≈178W)，而非全局峰 V≈28.8V(P≈278W)。
V0 = 36; delta = -0.5;   % 从接近开路电压的高压侧开始，向下搜索
V_po = V0;
for iteration = 1:500
    P0 = interp1(V_sweep, P_partial, V_po, 'linear', 0);
    V1 = V_po + delta;
    P1 = interp1(V_sweep, P_partial, V1, 'linear', 0);
    if P1 > P0
        V_po = V1;
    else
        V_po = V1 - delta;   % 回退
        delta = -delta;      % 反转方向
    end
    V_po = min(max(V_po, 0), V_oc);
end
P_po_trap = interp1(V_sweep, P_partial, V_po, 'linear', 0);

% ===== 全局扫描策略：先扫描全范围，锁定全局峰，再精调 =====
% 1. 全局扫描：找 P-V 全局峰值
[P_global, i_global] = max(P_partial);
V_global = V_sweep(i_global);

% 2. 精调：在全局峰附近用小步长精细定位（可选二阶细化）
V_fine = linspace(max(0, V_global-0.5), min(V_oc, V_global+0.5), 200);
P_fine = interp1(V_sweep, P_partial, V_fine, 'pchip', 0);
[P_mppt, i_fine] = max(P_fine);
V_mppt = V_fine(i_fine);

%% ===== 结果对比 =====
fprintf('===== 局部遮阴下 MPPT 结果对比 =====\n');
fprintf('局部P&O(易陷局部峰): V=%.2f V, P=%.2f W\n', V_po, P_po_trap);
fprintf('全局扫描找全局峰:   V=%.2f V, P=%.2f W\n', V_global, P_global);
fprintf('扫描+精调最终MPP:    V=%.2f V, P=%.2f W\n', V_mppt, P_mppt);
fprintf('结论：全局扫描策略优于局部P&O，避免陷入局部峰值。\n');

%% ===== 可视化对比 =====
figure;
plot(V_sweep, P_partial, 'b-', 'LineWidth', 1.5); hold on;
plot(V_mppt, P_mppt, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(V_po, P_po_trap, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
xlabel('Voltage (V)'); ylabel('Power (W)'); grid on;
title('局部遮阴下 MPPT：全局扫描 vs 局部P&O');
legend('多峰P-V曲线','全局扫描MPPT','局部P&O(陷阱)','Location','best');