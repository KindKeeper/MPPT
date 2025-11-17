function[V_PV,I_PV] = main_MP()
% 主函数：实现光伏电池最大功率点跟踪(MPPT)算法
% 输出：V_PV - 电压向量，I_PV - 电流向量

% ==== 光伏电池参数定义部分 ====
V_oc = 36;         % 开路电压（V），光伏电池在无负载条件下的最大电压
I_sc = 7.8;        % 短路电流（A），光伏电池在短路条件下的最大电流
I_o = 1e-7;        % 反向饱和电流[A]，二极管反向饱和电流，影响暗电流特性
q = 1.6e-19;       % 元电荷[C]，基本电荷常数
n = 1.5;           % 理想因子，反映二极管的非理想程度，典型范围1-2
k = 1.38e-23;      % 玻尔兹曼常数[J/K]，热力学常数
R_sh = 300;        % 并联电阻(Ω)，模拟电池的旁路漏电效应
T = 298;           % 温度（K），工作温度，对应25摄氏度
num_cells = 45;    % 光伏电池数量，串联的电池片数
P_threshold = 2;   % 功率变化阈值，用于判断是否远离最大功率点
V_threshold = 0.5; % 电压变化阈值，用于判断是否接近最大功率点

% 变步长系数定义
k1 = 0.02;         % 远离最大功率点时的步长增益系数
k2 = 0.005;        % 接近最大功率点时的步长增益系数

% ==== 初始化参数和数组 ====
% 生成从0到开路电压的1000个均匀分布的电压点
V_PV = linspace(0, V_oc, 1000);
% 预分配电流和功率数组
I_PV = zeros(size(V_PV));
P = zeros(size(V_PV));

% 光生电流计算，理想情况下等于短路电流
I_ph = I_sc;
    
% ==== 变步长扰动观察法主循环 ====
% 通过迭代调整工作点，跟踪最大功率点位置
for iter = 1:1000
    % 第一步初始化步长
    if iter == 1
       step = 0.001;  % 初始步长设置
    end
    
    % 计算二极管暗电流分量（肖克利二极管方程）
    I1 = I_o*(exp((q*V_PV)/(n*k*T*num_cells))-1);
    % 计算并联电阻漏电流
    I2 = (V_PV)/R_sh;
    % 总输出电流 = 光生电流 - 暗电流 - 漏电流
    I_PV = I_ph - I1 - I2;
    % 计算功率：P = V × I
    P = V_PV .* I_PV;
    
    % 预测下一个电压点
    V_new = V_PV + step;
    % 计算新电压点对应的电流值
    I_new = I_PV - I_o*(exp((q*V_new)/(n*k*T*num_cells))-1) - (V_new)/R_sh;
    % 计算新电压点对应的功率值
    P_new = V_new .* I_new;
    
    % 计算功率变化量和电压变化量的绝对值
    delta_P = abs(P_new - P);
    delta_V = abs(V_new - V_PV);    
    
    % ==== 变步长策略：根据工作点位置动态调整步长 ====
    if delta_P > P_threshold
        % 远离最大功率点区域，功率变化较大
        % 采用较大步长快速接近最大功率点
        step = k1 * delta_P;
    elseif delta_V < V_threshold
        % 接近最大功率点区域，电压变化较小
        % 采用较小步长精细调整，避免振荡
        % 根据功率变化方向决定步长符号（扰动观察法核心）
        step = sign(P_new - P) * k2 .* delta_P;
    else
        % 中间状态，使用默认步长
        step = 0.001;
    end
    
    % 更新电压向量为新的预测值
    V_PV = V_new;

    % ==== 收敛判断：功率变化小于阈值时终止迭代 ====
    if abs(delta_P) < 1e-3
        % 功率变化足够小，认为已收敛到最大功率点
        break;
    end
    
    % 记录当前迭代的最大功率点电压和功率
    V_mppt = max(V_PV);
    P_mppt = max(P);
end

% 确保电流值为非负（物理约束）
I_PV = max(0, I_PV);

% ==== 结果可视化部分 ====

% 绘制I-V特性曲线（电流-电压关系）
figure;
plot(V_PV, I_PV, '-b');  % 蓝色实线
xlabel('Voltage (V)')    % X轴标签：电压
ylabel('Current (A)')    % Y轴标签：电流

% 绘制P-V特性曲线（功率-电压关系），红色星号标记
figure;
plot(V_PV, P, '*r');    % 红色星号
xlabel('Voltage (V)')    % X轴标签：电压
ylabel('Power (W)')      % Y轴标签：功率

% 在命令窗口显示最大功率点结果
disp(['最大功率点电压：',num2str(V_mppt),'V,最大功率:',num2str(P_mppt),'W'])

% 再次绘制P-V曲线，使用红色虚线
figure;
plot(V_PV, P, '--r');   % 红色虚线
xlabel('V (V)')          % X轴标签：电压
ylabel('p (w)')          % Y轴标签：功率

end
