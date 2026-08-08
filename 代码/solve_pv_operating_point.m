function [Vsol, Isol] = solve_pv_operating_point(D, Rload, G, T, params)
% 功能：求解给定占空比下的光伏工作点电压和电流
% 输入：D-占空比, Rload-负载电阻, G-辐照度, T-温度, params-光伏参数
% 输出：Vsol-工作点电压, Isol-工作点电流

% 计算Boost变换器在光伏侧看到的等效输入电阻
% 理想Boost变换器关系：Rin = Rload * (1-D)^2
Rin = Rload * (1 - D)^2;
if Rin <= 0, Rin = 1e-6; end  % 防零除保护

% 定义残差方程：光伏输出电流等于负载电流
% f(V) = Ipv(V) - V/Rin = 0
f = @(V) pv_current_given_voltage(max(V,0), G, T, params) - max(V,0)/Rin;

% 估算开路电压作为搜索区间上限（简化温度修正）
Voc_est = params.Voc_stc * (1 + 0.0*(T-25)/25); % 简化处理，实际应加温度系数

% TODO: 开路电压温度修正
% 问题：未考虑温度对开路电压的影响，实际Voc随温度升高而降低
% 建议：添加温度系数，如 Voc_est = params.Voc_stc + beta*(T-25)，beta ≈ -0.003/°C

V0 = 0.7*Voc_est;  % 初始猜测值，设为估计开路电压的70%

% TODO: fzero初始猜测优化
% 问题：固定比例初始值在极端条件下（如低辐照度）可能不理想
% 建议：基于当前工作条件动态调整初始猜测值

% 使用fzero求解非线性方程
try
    % 优先使用区间搜索，提高收敛可靠性
    Vsol = fzero(f, [0, max(1e-6, 1.2*Voc_est)]);
catch
    % 区间搜索失败时使用单点搜索
    Vsol = fzero(f, V0);
end
Vsol = max(Vsol, 0);  % 确保电压非负
Isol = pv_current_given_voltage(Vsol, G, T, params);  % 计算对应电流
