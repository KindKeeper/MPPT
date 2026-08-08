function [Vvec, Ivec] = sweep_pv_iv(G, T, params)
% 功能：扫描生成完整I-V特性曲线
% 输入：G-辐照度, T-温度, params-光伏参数
% 输出：Vvec-电压向量, Ivec-电流向量

Vmax = params.Voc_stc*1.15;  % 设置扫描电压上限（略大于开路电压）
Vvec = linspace(0, Vmax, 200);  % 生成200个电压点
Ivec = zeros(size(Vvec));       % 预分配电流数组

% 遍历所有电压点计算对应电流
for i = 1:numel(Vvec)
    Ivec(i) = pv_current_given_voltage(Vvec(i), G, T, params);
end
