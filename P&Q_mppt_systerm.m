function[V_PV,I_PV] = main_MP()
V_oc = 36; % 开路电压（V）
I_sc = 7.8; % 短路电流（A）
I_o = 1e-7;   % 反向饱和电流 [A]
q = 1.6e-19;  % 元电荷 [C]
n = 1.5; % 理想因子
k = 1.38e-23; % 玻尔兹曼常数 [J/K]
R_sh = 300; % 并联电阻(Ω)
T = 298; % 温度（K）
num_cells = 45; % 光伏电池数量
P_threshold = 2; % 设定阈值
V_threshold = 0.5; % 设定阈值
k1 = 0.02;
k2 = 0.005;


% 参数范围设置
V_PV = linspace(0, V_oc, 1000);
I_PV = zeros(size(V_PV));
P = zeros(size(V_PV));

% 光生电流
I_ph = I_sc;
    
% 输出电流与输出电压的计算   
for iter = 1:1000
    if iter == 1
       step = 0.001;
    end
    I1 = I_o*(exp((q*V_PV)/(n*k*T*num_cells))-1);
    I2 = (V_PV)/R_sh;
    I_PV = I_ph - I1 - I2;
    P = V_PV .* I_PV;
    V_new = V_PV + step;
    I_new = I_PV-I_o*(exp((q*V_new)/(n*k*T*num_cells))-1)-(V_new)/R_sh;
    P_new = V_new .* I_new;
    delta_P = abs(P_new - P);
    delta_V = abs(V_new - V_PV);    
    
% 判断功率变化方向
if delta_P > P_threshold
   step = k1 * delta_P ; % 大功率下步长变大 
   elseif delta_V < V_threshold
          step = sign(P_new - P) * k2 .* delta_P; % 接近极值点时缩小步长
   else
   step = 0.001; % 默认步长
end
% 更新电压
V_PV = V_new;

% 终止条件: 功率变化小于阈值
if abs(delta_P) < 1e-3
   break;
end
V_mppt = max(V_PV);
P_mppt = max(P);
end
I_PV = max(0, I_PV); % 确保电流非负
figure;
plot(V_PV, I_PV, '-b');
xlabel('Voltage (V)')
ylabel('Current (A)')
figure;
plot(V_PV, P, '*r');
xlabel('Voltage (V)')
ylabel('Power (W)')
disp(['最大功率点电压：',num2str(V_mppt),'V,最大功率:',num2str(P_mppt),'W'])   
% 可视化结果
figure;
plot(V_PV, P, '--r');
xlabel('V (V)')
ylabel('p (w)')
end
