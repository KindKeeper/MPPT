% 扰动观察法I_PV
% 光伏电池参数
I_o = 1e-7;   % 反向饱和电流 [A]
q = 1.6e-19;  % 元电荷 [C]
T = 298;      % 室温 [K]
n = 1.5;      % 理想因子
k = 1.38e-23; % 玻尔兹曼常数 [J/K]
step = 0.01; % 扰动步长 [V]
k1 = 0.03;
k2 = 0.005;
P_threshold = 2; % 设定阈值
V_threshold = 0.5; % 设定阈值
% P&O算法主循环
for iter = 1:10000

    % 计算当前功率
    I = I_PV - I_o * (exp((q * V_PV)/(n * k * T)) - 1);
    P = V_PV * I;

    % 电压扰动
    V_new = V_PV + step;
    I_new = I_PV - I_o * (exp((q * V_new)/(n * k * T)) - 1);
    P_new = V_new * I_new;
    delta_P = abs(P_new - P);
    delta_V = abs(V_new - V_PV);    
    
    % 判断功率变化方向
    if delta_P > P_threshold
        step = k1 * delta_P ; % 大功率下步长变大 
    elseif delta_V < V_threshold
            step = -k2 * delta_P; % 接近极值点时缩小步长
    else
        step = step; % 默认步长
    end

    % 更新电压
    V_PV = V_new;

    % 终止条件: 功率变化小于阈值
    if abs(P_new - P) < 1e-3
       break;
    end
end

V_mppt = V_PV;
P_mppt = P;