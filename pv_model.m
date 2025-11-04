function[V_PV,I_PV] = PV_model()
V_oc = 36.7; % 开路电压（V）
I_sc = 6.7; % 短路电流（A）
I_o = 3e-10;   % 反向饱和电流 [A]
q = 1.6e-19;  % 元电荷 [C]
n = 1; % 理想因子
k = 1.38e-23; % 玻尔兹曼常数 [J/K]
R_sh = 300; % 并联电阻(Ω)
T = 298; % 温度（K）
num_cells = 60; % 光伏电池数量

% 参数范围设置
V_PV = linspace(0, V_oc, 1000);
I_PV = zeros(size(V_PV));
P = zeros(size(V_PV));

% 光生电流
I_ph = I_sc;
    
% 输出电流与输出电压的计算   
for iter = 1:1000
    I1 = I_o*(exp((q*V_PV)/(n*k*T*num_cells))-1);
    I2 = V_PV / R_sh;
    I_PV = I_ph - I1 -I2;
    P = V_PV .* I_PV;    
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
end