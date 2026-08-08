function I = pv_current(V, pv, G, T)
% pv_current.m — 光伏单二极管模型电流计算（全仓库统一模型）
% 功能：给定电压，用单二极管模型计算光伏组件输出电流
% 输入：
%   V — 电压 [V]（标量或向量）
%   pv — 参数结构体（由 pv_params() 生成）
%   G — 辐照度 [W/m^2]，默认 1000（STC）
%   T — 组件温度 [°C]，默认 25（STC）
% 输出：
%   I — 电流 [A]（与 V 同形状，物理上截断为非负）
% 运行环境：MATLAB 或 Octave
% 说明：全仓库算法的 PV 模型统一入口，保证各脚本结果可比。

if nargin < 3, G = 1000; end
if nargin < 4, T = 25;   end

V = max(V, 0);

% 热电压（考虑串联片数）
Vt = pv.Ns * pv.k * (T + 273.15) / pv.q;

% 光生电流：随辐照度线性变化 + 温度微调
Iph = pv.Isc_stc * (G/1000) + pv.alpha_Isc * (T - 25);

% 反向饱和电流（由 STC 条件反推）
Vt_stc = pv.Ns * pv.k * (25 + 273.15) / pv.q;
I0_stc = pv.Isc_stc / (exp(pv.Voc_stc/(pv.n*Vt_stc)) - 1);

% I0 温度修正（含禁带宽度 Eg 影响）
T_K = T + 273.15;
I0 = I0_stc * (T_K/298.15).^3 .* ...
    exp((pv.q*pv.Eg/pv.k) * (1/298.15 - 1/T_K) / (pv.n*pv.Ns));

% 单二极管隐式方程解（忽略并联电阻的简化解析式，含 Rs 项近似）
% I = Iph - I0*(exp((V + I*Rs)/(n*Vt)) - 1) - (V + I*Rs)/Rsh
% 采用固定点迭代求解 I（对典型组件收敛快速、数值稳定）
I = Iph;
for iter = 1:30
    I = Iph - I0*(exp((V + I.*pv.Rs)./(pv.n*Vt)) - 1) - (V + I.*pv.Rs)./pv.Rsh;
    I = max(I, 0);
end

end