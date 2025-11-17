%% 光伏参数
V_oc = 36.6;       % 开路电压(V)
I_sc = 6.2;        % 短路电流(A)
I_o  = 3e-10;      % 反向饱和电流 (A)
q = 1.6e-19;       % 元电荷 (C)
n = 1;             % 理想因子
k = 1.38e-23;      % 玻尔兹曼常数 (J/K)
R_sh = 300;        % 并联电阻(Ω)
T = 298;           % 温度(K)
num_c = 60;        % 串联电池数目
P_threshold = 2;   % 功率阈值
V_threshold = 1;   % 电压阈值

% 变步长系数
k1 = 0.2;          % 远离峰值
k2 = 0.02;         % 接近峰值

V_PV = linspace(0, V_oc, 1000).';
I_ph = I_sc;
I1 = I_o * (exp((q*V_PV)./(n*k*T*num_c)) - 1);
I2 = V_PV / R_sh;
I_PV = I_ph - I1 - I2;
I_PV = max(I_PV, 0);  % 保证非负
P = V_PV .* I_PV;

S = table(V_PV, I_PV, P, 'VariableNames', {'V','I','P'}); % S 为数据表格

N   = height(S);
dV  = mean(diff(S.V)); % 把步长反应为格数
iter = 1000;
idx = 1;

% 记录迭代过程的数据
V_s = zeros(1,iter);  % 电压数据
I_s = zeros(1,iter);  % 电流数据
P_s = zeros(1,iter);  % 功率数据
for k = 1:iter

    % 更新电压
    V = S.V(idx);   P = S.P(idx);
    V_new = S.V(idx + 1); P_new = S.P(idx + 1);
    dP = P_new - P;
    dV = V_new - V; 

    % 记录电压
    V_s(k) = V;P_s(k) = P;

    % 计算斜率 dP/dV 的符号
    sgn = sign(dP * dV);

    % 确认步长
    if abs(dP) > P_threshold
        step = k1 * abs(dP);  % 远离峰值，大步
    elseif abs(dV) < V_threshold
        step = k2 * abs(dP);  % 接近峰值，小步
    else
        step = 0.1;  % 默认步长
    end

    % 将电压步长转化为表格的跨幅
    step_idx = max(1, round(step / max(dV, 1e-9)));
    idx_new = idx + sgn * step_idx;
    idx = min(max(idx_new, 1), N-1); % 确定下一次取哪一格

end

V_mppt = S.V(idx);
P_mppt = S.P(idx);

figure; 
plot(S.V, S.I, '-b');
xlabel('Voltage (V)'); 
ylabel('Current (A)');
title('I–V'); 
hold on; plot(V_mppt, S.I(idx), 'ro');

figure; 
plot(S.V, S.P, '-r');
xlabel('Voltage (V)'); 
ylabel('Power (W)');
title('P–V'); 
hold on; plot(V_mppt, P_mppt, 'ko');

figure; 
plot(P_s,'LineWidth',1.2);
xlabel('迭代过程'); 
ylabel('P (W)');

figure; 
plot(V_s,'LineWidth',1.2);  
xlabel('迭代过程'); 
ylabel('V_{pv} (V)');


 
