% 扰动观察法
% 光伏电池参数
Iph = 10;     % 光生电流 [A]
I_o = 1e-7;   % 反向饱和电流 [A]
q = 1.6e-19;  % 元电荷 [C]
T = 298;      % 室温 [K]
n = 1.5;      % 理想因子
k = 1.38e-23; % 玻尔兹曼常数 [J/K]
delta = 0.01; % 扰动步长 [V]
V = 10;       % 初始电压 [V]

% P&O算法主循环
for iter = 1:10000

    % 计算当前功率
    I = Iph - I_o*(exp((q*V)/(n*k*T)) - 1);
    P = V * I;

    % 电压扰动
    V_new = V + delta;
    I_new = Iph - I_o*(exp((q*V_new)/(n*k*T)) - 1);
    P_new = V_new*I_new;
    dP = P_new - P;
    dV = V_new - V;    
    
    % 判断功率变化方向
    if P_new > P
        delta = delta * abs(dP/dV); % 变步长正向扰动 
    else
        delta = -delta * abs(dP/dV); % 变步长反向扰动
    end

    % 更新电压
    V = V_new;

    % 终止条件: 功率变化小于阈值
    if abs(P_new - P) < 1e-3
       break;
    end
end

% 输出MPPT结果
disp(['MPPT电压: ', num2str(V), 'V, 功率: ', num2str(P), 'W']);