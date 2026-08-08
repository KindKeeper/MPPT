function [V_PV, I_PV] = pv_single_diode_model()
% 光伏电池单二极管模型函数
% 功能：计算并绘制光伏电池的电流-电压(I-V)和功率-电压(P-V)特性曲线
% 输出：V_PV - 电压向量，I_PV - 电流向量
% 运行环境：MATLAB 或 Octave
% 用法：直接调用 pv_single_diode_model 即可看到曲线和最大功率结果。
%       光伏模型参数统一来自 pv_params() + pv_current()。

% ==== 光伏组件参数（统一参数源） ====
pv = pv_params();
V_oc = pv.Voc_stc;

% ==== 生成电压采样点并计算 I-V/P-V 特性 ====
V_PV = linspace(0, V_oc, 1000);
I_PV = pv_current(V_PV, pv);   % 统一单二极管模型
P = V_PV .* I_PV;               % 输出功率

% 最大功率点
[P_mppt, imax] = max(P);
V_mppt = V_PV(imax);

fprintf('开路电压: %.2f V, 短路电流: %.2f A\n', V_oc, pv.Isc_stc);
fprintf('最大功率点: V=%.2f V, P=%.2f W\n', V_mppt, P_mppt);

% ==== 结果可视化 ====
% I-V 特性曲线
figure;
plot(V_PV, I_PV, '-b');
xlabel('Voltage (V)');
ylabel('Current (A)');
title('I-V 特性曲线');
grid on;

% P-V 特性曲线
figure;
plot(V_PV, P, '-r');
xlabel('Voltage (V)');
ylabel('Power (W)');
title('P-V 特性曲线');
grid on;
hold on;
plot(V_mppt, P_mppt, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');  % 标记MPP

end