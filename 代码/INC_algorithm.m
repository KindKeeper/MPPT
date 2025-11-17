% 电导增量法
% 光伏电池参数
Voc = 36; % 开路电压
Isc = 5;  % 短路电流
Vmpp = 35; % 最大功率点电压
Pmax = 800; % 最大功率
delta = 0.01; % 扰动步长
tolerance = 1e-6; % 收敛容差

Vpv = Voc/2; % 初始电压猜测
while true
    Ipv = Isc - (Isc/Voc)*Vpv; % 简化电流模型
    dI_dV = -Isc/Voc; % 电流对电压导数
    dP_dV = Ipv + Vpv*dI_dV; % 功率对电压导数
    if abs(dP_dV) > tolerance
        if dP_dV > 0            
            Vpv += delta; % 沿导纳增加方向调整
        else
            Vpv -= delta;
        end
    else
        break; % 收敛到MPP
    end
end

disp(['最大功率点电压: ', num2str(Vpv)]);