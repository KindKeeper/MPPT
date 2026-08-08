function pv = pv_params(varargin)
% pv_params.m — 光伏组件标准参数（统一参数源）
% 功能：返回一套标准 PV 组件的物理参数结构体，供全仓库所有算法脚本共用。
%       统一参数是"可复现"的第一步：一处修改，全仓库生效。
% 用法：pv = pv_params();                 % 默认标准组件
%       pv = pv_params('组件名');          % 预留：未来支持多型号
% 输出结构体字段：
%   pv.Ns        — 串联电池片数
%   pv.Isc_stc   — STC 短路电流 [A]
%   pv.Voc_stc   — STC 开路电压 [V]
%   pv.n         — 二极管理想因子
%   pv.Rs        — 串联电阻 [Ohm]
%   pv.Rsh       — 并联电阻 [Ohm]
%   pv.alpha_Isc — 电流温度系数 [A/°C]
%   pv.Eg        — 硅禁带宽度 [eV]
%   pv.q, pv.k   — 元电荷、玻尔兹曼常数（SI）
% 运行环境：MATLAB 或 Octave

% 默认：常见 60 片单晶硅组件（与 mppt_adaptive_slope_boost 保持一致的量级）
pv.Ns        = 60;
pv.Isc_stc   = 8.7;    % 短路电流 [A]
pv.Voc_stc   = 37.0;   % 开路电压 [V]
pv.n         = 1.3;    % 理想因子
pv.Rs        = 0.35;   % 串联电阻 [Ohm]
pv.Rsh       = 500;    % 并联电阻 [Ohm]
pv.alpha_Isc = 0.005;  % 电流温度系数 [A/°C]
pv.Eg        = 1.12;   % 禁带宽度 [eV]

% 物理常数
pv.q = 1.602176634e-19;   % 元电荷 [C]
pv.k = 1.380649e-23;      % 玻尔兹曼常数 [J/K]

end