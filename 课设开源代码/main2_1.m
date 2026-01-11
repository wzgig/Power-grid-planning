% 最终推荐方案的确定程序一
clc,clear,close all
tic
format long g

%% =========================
%  节点数量及距离矩阵初始化
%% =========================
number_of_point = 5;
d = zeros(number_of_point);
d(1,2) = 130; d(1,3) = 260; d(1,4) = 160; d(1,5) = 160;
d(2,3) = 130; d(2,4) = 140; d(2,5) = 140;
d(3,4) = 160; d(3,5) = 160;
d(4,5) = 280;
d = d + d.';   % 对称化

%% =========================
%  基础参数（仅保留 A、负荷、Slack）
%% =========================
S_A = 680 + 1j*510;          % 发电厂A两台机组满载（视为"总额定"）
k3_A = 10/100;               % 厂用电率（只扣P）

S_S = 0 + 0j;
S1_min = 95 + 1j*71.25;
S2_min = 85 + 1j*41.17;
S3_min = 90 + 1j*55.77;
SA_min = 24 + 1j*18;
S1_max = 270 + 1j*202.5;
S2_max = 260 + 1j*125.92;
S3_max = 280 + 1j*173.52;
SA_max = 70  + 1j*52.5;      % A点直配负荷（注意：这是负荷，不是发电出力）

% ====== 负荷（正值表示负荷）======
Pd1 = abs(real(S1_max));  Qd1 = abs(imag(S1_max));
Pd2 = abs(real(S2_max));  Qd2 = abs(imag(S2_max));
Pd3 = abs(real(S3_max));  Qd3 = abs(imag(S3_max));
PdA_load = abs(real(SA_max));  QdA_load = abs(imag(SA_max));

% ====== Slack（无穷大系统）负荷固定为0 ======
PdS = 0; QdS = 0;

%% =========================
%  阻抗及基准值
%% =========================
z1 = 0.05 + 1j*0.4;
Z_all = z1 .* d;

SB = 100;
UB = 220;
ZB = UB.^2 ./ SB;

% 清理无用变量
clear S1_min S2_min S3_min S1_max S2_max S3_max SA_max

%%% 方案 5 计算
disp('方案 5 潮流:')
SA_case_list = [S_A; 0.5*S_A];
qk = 1;
SA_case = SA_case_list(qk);
% A发电机出力（P扣厂用电，Q不扣）
PdA = real(SA_case) * (1 - k3_A);  % 这里沿用你 Plan_1 参数名 PdA/QdA
QdA = imag(SA_case);

[U5,theta5,P5,Q5,Branch_Loss5] = Plan_5( ...
    Pd1, Pd2, Pd3,  Qd1, Qd2, Qd3, ...
    Z_all, ZB, UB, SB, ...
    PdA, QdA, PdA_load, QdA_load, ...
    number_of_point);
SIJ_5 = Branch_Loss5(:,end-1) + 1j*Branch_Loss5(:,end);
S5 = 1000*(abs(SIJ_5))./sqrt(3)/220/1.15;


%%% 方案 10 计算
disp('方案 10 潮流:')
SA_case_list = [S_A; 0.5*S_A];
qk = 1;
SA_case = SA_case_list(qk);

% A发电机出力（P扣厂用电，Q不扣）
PdA = real(SA_case) * (1 - k3_A);  % 这里沿用你 Plan_1 参数名 PdA/QdA
QdA = imag(SA_case);

[U10,theta10,P10,Q10,Branch_Loss10] = Plan_10( ...
    Pd1, Pd2, Pd3,  Qd1, Qd2, Qd3, ...
    Z_all, ZB, UB, SB, ...
    PdA, QdA, PdA_load, QdA_load, ...
    number_of_point);
SIJ_10 = Branch_Loss10(:,end-1) + 1j*Branch_Loss10(:,end);
S10 = 1000*(abs(SIJ_10))./sqrt(3)/220/1.15;

% 导线临界电压计算
[~,~,~,U_cr1] = calculate_daoxian('LGJ-240',1);
Uphi1 = 110/sqrt(3);
[~,~,~,U_cr2] = calculate_daoxian('LGJ-400',1);
Uphi2 = 220/sqrt(3);

% 方案 5 导线阻抗导纳计算
R_Plan5 = zeros(number_of_point);
X_Plan5 = zeros(number_of_point);
B_Plan5 = zeros(number_of_point);
[R_Plan5(1,4),X_Plan5(1,4),B_Plan5(1,4),~] = ...
    calculate_daoxian('LGJ-400',d(1,4));
[R_Plan5(1,5),X_Plan5(1,5),B_Plan5(1,5),~] = ...
    calculate_daoxian('LGJ-240',d(1,5));
[R_Plan5(2,4),X_Plan5(2,4),B_Plan5(2,4),~] = ...
    calculate_daoxian('LGJ-400',d(2,4));
[R_Plan5(2,5),X_Plan5(2,5),B_Plan5(2,5),~] = ...
    calculate_daoxian('LGJ-400',d(2,5));
[R_Plan5(3,4),X_Plan5(3,4),B_Plan5(3,4),~] = ...
    calculate_daoxian('LGJ-400',d(3,4));
[R_Plan5(3,5),X_Plan5(3,5),B_Plan5(3,5),~] = ...
    calculate_daoxian('LGJ-400',d(3,5));
R_Plan5 = R_Plan5 + R_Plan5.';
X_Plan5 = X_Plan5 + X_Plan5.';
B_Plan5 = B_Plan5 + B_Plan5.';
B_Plan5(B_Plan5==0) = 1e5;
Daoxian_Plan5 = [
    R_Plan5(1,4),X_Plan5(1,4),B_Plan5(1,4);...
    R_Plan5(1,5),X_Plan5(1,5),B_Plan5(1,5);...
    R_Plan5(2,4),X_Plan5(2,4),B_Plan5(2,4);...
    R_Plan5(2,5),X_Plan5(2,5),B_Plan5(2,5);...
    R_Plan5(3,4),X_Plan5(3,4),B_Plan5(3,4);...
    R_Plan5(3,5),X_Plan5(3,5),B_Plan5(3,5)
];

% 方案 10 导线阻抗导纳计算
R_Plan10 = zeros(number_of_point);
X_Plan10 = zeros(number_of_point);
B_Plan10 = zeros(number_of_point);
% 节点对(1,2) - 导线型号LGJ-240（沿用原有规格）
[R_Plan10(1,2),X_Plan10(1,2),B_Plan10(1,2),~] = calculate_daoxian('LGJ-240',d(1,2));

% 节点对(1,4) - 导线型号LGJ-400（新增）
[R_Plan10(1,4),X_Plan10(1,4),B_Plan10(1,4),~] = calculate_daoxian('LGJ-400',d(1,4));

% 节点对(1,5) - 导线型号LGJ-400（保留原有）
[R_Plan10(1,5),X_Plan10(1,5),B_Plan10(1,5),~] = calculate_daoxian('LGJ-400',d(1,5));

% 节点对(2,3) - 导线型号LGJ-400（新增）
[R_Plan10(2,3),X_Plan10(2,3),B_Plan10(2,3),~] = calculate_daoxian('LGJ-400',d(2,3));

% 节点对(2,4) - 导线型号LGJ-400（保留原有）
[R_Plan10(2,4),X_Plan10(2,4),B_Plan10(2,4),~] = calculate_daoxian('LGJ-400',d(2,4));

% 节点对(2,5) - 导线型号LGJ-400（保留原有）
[R_Plan10(2,5),X_Plan10(2,5),B_Plan10(2,5),~] = calculate_daoxian('LGJ-400',d(2,5));

% 节点对(3,4) - 导线型号LGJ-400（新增）
[R_Plan10(3,4),X_Plan10(3,4),B_Plan10(3,4),~] = calculate_daoxian('LGJ-400',d(3,4));
R_Plan10 = R_Plan10 + R_Plan10.';
X_Plan10 = X_Plan10 + X_Plan10.';
B_Plan10 = B_Plan10 + B_Plan10.';
B_Plan10(B_Plan10==0) = 1e5;
Daoxian_Plan10 = [
    R_Plan10(1,2),X_Plan10(1,2),B_Plan10(1,2);...
    R_Plan10(1,4),X_Plan10(1,4),B_Plan10(1,4);...
    R_Plan10(1,5),X_Plan10(1,5),B_Plan10(1,5);...
    R_Plan10(2,3),X_Plan10(2,3),B_Plan10(2,3);...
    R_Plan10(2,4),X_Plan10(2,4),B_Plan10(2,4);...
    R_Plan10(2,5),X_Plan10(2,5),B_Plan10(2,5);...
    R_Plan10(3,4),X_Plan10(3,4),B_Plan10(3,4)
];

% 整合并写入输电线参数文件
SDX = [R_Plan5,X_Plan5,B_Plan5;R_Plan10,X_Plan10,B_Plan10];
writematrix(SDX,'输电线参数.txt');

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数：计算导线的电阻、电抗、电纳和临界电压
% 输入：type-导线型号(LGJ-240/LGJ-400)，l-长度(km)
% 输出：R-电阻(Ω)，X-电抗(Ω)，B-电纳(S)，U_cr-临界电压(kV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R,X,B,U_cr] = calculate_daoxian(type,l)
    if strcmp(type,'LGJ-240')
        C_L = 0.905*1e-8;
        S = 240;
    elseif strcmp(type,'LGJ-400')
        C_L = 0.87*1e-8;
        S = 400;
    end
    rho = 31.5;
    r = sqrt(S/pi);
    lgDmcr = 7.56*1e-6/100/pi/C_L;
    R = l*rho/S;
    X = l*100*pi*(4.6*lgDmcr+0.5)*1e-4;
    B = l*7.56*1e-6/lgDmcr;
    U_cr = 49.3*0.9*r/10*lgDmcr;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数：潮流计算核心函数
% 输入：mpc-潮流计算结构体
% 输出：U-节点电压幅值，theta-节点电压相角，P-节点注入有功，Q-节点注入无功，Branch_Loss-支路损耗
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U,theta,P,Q,Branch_Loss] = Power_Flow_Calculation(mpc)
    Power_flow_result = runpf(mpc,mpoption('pf.nr.max_it',100));
    U = Power_flow_result.bus(:,8);
    theta = Power_flow_result.bus(:,9);
    
    % 节点注入功率-->输入节点为+,输出为负
    AA = zeros(size(Power_flow_result.bus,1),1);
    BB = AA;
    idx = (Power_flow_result.bus(:,2)~=1);
    number = Power_flow_result.bus(idx,1);
    AA(number) = Power_flow_result.gen(:,2);
    BB(number) = Power_flow_result.gen(:,3);
    P = AA - Power_flow_result.bus(:,3);
    Q = BB - Power_flow_result.bus(:,4);
    
    % 支路损耗： 起点,终点,有功网损,无功网损,输送功率,输送无功功率
    Branch_Loss = [
        Power_flow_result.branch(:,1),Power_flow_result.branch(:,2),...
        Power_flow_result.branch(:,14)+Power_flow_result.branch(:,16),...
        Power_flow_result.branch(:,15)+Power_flow_result.branch(:,17),...
        Power_flow_result.branch(:,14),Power_flow_result.branch(:,15)
    ];
    % 以传入线路的有功功率绝对值作为输送功率
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 方案5计算函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U,theta,P,Q,Branch_Loss] = Plan_5( ...
    Pd1, Pd2, Pd3,  ...
    Qd1, Qd2, Qd3,  ...
    Z_all, ZB, UB, SB, ...
    PdA, QdA, ...
    PdA_load, QdA_load, ...
    number_of_point)

f = zeros(number_of_point); % 节点数×节点数的全0矩阵
a = ones(number_of_point);% 全1矩阵，节点之间是否双回路接线
% 定义支路连接：f(i,j)=1表示节点i-j有支路
f(1,4) = 1; f(1,5) = 1;
f(2,4) = 1; f(2,5) = 1;
f(3,4) = 1; f(3,5) = 1;
a(1,4) = 2;
a(2,4) = 2;
f = f + f.';
a = a + a.';
% 计算标幺值阻抗矩阵
z = Z_all./a./f./ZB;
z(1:length(z)+1:end) = 0;
z(real(z)==inf) = 1e5*(1+1j);
r = real(z);
x = imag(z);

% 构建潮流计算MPC结构体
mpc = struct();
mpc.version = '2';
mpc.baseMVA = SB; % 功率/容量基准值

%% 母线数据输入
% 节点编号 类型 Pd Qd Gs Bs area Vm Vang baseKV zone Vmax Vmin
mpc.bus = [
    1 1 Pd1 Qd1 0 0 1 1 0 UB 1 1.1 0.94;
    2 1 Pd2 Qd2 0 0 1 1 0 UB 1 1.1 0.94;
    3 1 Pd3 Qd3 0 0 1 1 0 UB 1 1.1 0.94;
    4 1 -(PdA-PdA_load) -(QdA-QdA_load) 0 0 1 1 0 UB 1 1.1 0.94;
    5 3 0 0 0 0 1 1 0 UB 1 1.1 0.94;
    ];

%% 发电机组数据输入（无穷大电源专属配置）
% 节点编号 Pg Qg Qmax Qmin Vg mBase status Pmax Pmin
mpc.gen = [
    5 0 0 1e9 -1e9 1.0 SB 1 1e9 -1e9;  % 无穷大电源（核心修改）
    ];

%% 支路数据输入
% 起点 终点 r x b rateA rateB rateC 变比 angle status angmin angmax
mpc.branch = [
    % 第一组：f矩阵定义的6条有效支路，按节点号排序
    1 4 r(1,4) x(1,4) 0 0 0 0 0 0 1 -360 360;  % 支路1-2（匹配f(1,2)=1）
    1 5 r(1,5) x(1,5) 0 0 0 0 0 0 1 -360 360;  
    2 4 r(2,4) x(2,4) 0 0 0 0 0 0 1 -360 360;  % 支路2-4（匹配f(2,4)=1）
    2 5 r(2,5) x(2,5) 0 0 0 0 0 0 1 -360 360;  % 支路2-5（匹配f(2,5)=1）
    3 4 r(3,4) x(3,4) 0 0 0 0 0 0 1 -360 360;  % 支路3-4（匹配f(3,4)=1）
    3 5 r(3,5) x(3,5) 0 0 0 0 0 0 1 -360 360;
    ];

% 潮流计算与参数计算
[U,theta,P,Q,Branch_Loss] = Power_Flow_Calculation(mpc);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 方案10计算函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U,theta,P,Q,Branch_Loss] = Plan_10( ...
    Pd1, Pd2, Pd3,  ...
    Qd1, Qd2, Qd3,  ...
    Z_all, ZB, UB, SB, ...
    PdA, QdA, ...
    PdA_load, QdA_load, ...
    number_of_point)

f = zeros(number_of_point); % 节点数×节点数的全0矩阵
a = ones(number_of_point);% 全1矩阵，节点之间是否双回路接线
% 定义支路连接：f(i,j)=1表示节点i-j有支路
f(1,2) = 1; f(1,4) = 1; f(1,5) = 1;
f(2,4) = 1; f(2,5) = 1;
f(2,3) = 1;
f(3,4) = 1;
a(1,4) = 2; a(1,5) = 2;
a(2,4) = 1; a(2,5) = 1;
f = f + f.';
a = a + a.';
% 计算标幺值阻抗矩阵
z = Z_all./a./f./ZB;
z(1:length(z)+1:end) = 0;
z(real(z)==inf) = 1e5*(1+1j);
r = real(z);
x = imag(z);

% 构建潮流计算MPC结构体
mpc = struct();
mpc.version = '2';
mpc.baseMVA = SB; % 功率/容量基准值

%% 母线数据输入
% 节点编号 类型 Pd Qd Gs Bs area Vm Vang baseKV zone Vmax Vmin
mpc.bus = [
    1 1 Pd1 Qd1 0 0 1 1 0 UB 1 1.1 0.94;
    2 1 Pd2 Qd2 0 0 1 1 0 UB 1 1.1 0.94;
    3 1 Pd3 Qd3 0 0 1 1 0 UB 1 1.1 0.94;
    4 1 -(PdA-PdA_load) -(QdA-QdA_load) 0 0 1 1 0 UB 1 1.1 0.94;
    5 3 0 0 0 0 1 1 0 UB 1 1.1 0.94;
    ];

%% 发电机组数据输入（无穷大电源专属配置）
% 节点编号 Pg Qg Qmax Qmin Vg mBase status Pmax Pmin
mpc.gen = [
    5 0 0 1e9 -1e9 1.0 SB 1 1e9 -1e9;  % 无穷大电源（核心修改）
    ];

%% 支路数据输入
% 起点 终点 r x b rateA rateB rateC 变比 angle status angmin angmax
mpc.branch = [
    % 第一组：f矩阵定义的6条有效支路，按节点号排序
    1 2 r(1,2) x(1,2) 0 0 0 0 0 0 1 -360 360;  % 支路1-2（匹配f(1,2)=1）
    1 4 r(1,4) x(1,4) 0 0 0 0 0 0 1 -360 360;  
    1 5 r(1,5) x(1,5) 0 0 0 0 0 0 1 -360 360;  
    2 3 r(2,3) x(2,3) 0 0 0 0 0 0 1 -360 360;
    2 4 r(2,4) x(2,4) 0 0 0 0 0 0 1 -360 360;  % 支路2-4（匹配f(2,4)=1）
    2 5 r(2,5) x(2,5) 0 0 0 0 0 0 1 -360 360;  % 支路2-5（匹配f(2,5)=1）
    3 4 r(3,4) x(3,4) 0 0 0 0 0 0 1 -360 360;
    ];

% 潮流计算与参数计算
[U,theta,P,Q,Branch_Loss] = Power_Flow_Calculation(mpc);
end