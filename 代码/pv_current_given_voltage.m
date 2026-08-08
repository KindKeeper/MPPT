function I = pv_current_given_voltage(V, G, T, params)
% pv_current_given_voltage — 单二极管模型电流计算（薄封装）
% 功能：计算给定电压下的光伏输出电流（兼容 mppt_adaptive_slope_boost 接口）
% 输入：V-电压, G-辐照度, T-温度, params-光伏参数
% 输出：I-电流
% 说明：统一模型的封装，内部调用 pv_current.m，全仓库共用一套 PV 模型。
% 运行环境：MATLAB 或 Octave
I = pv_current(V, params, G, T);
end