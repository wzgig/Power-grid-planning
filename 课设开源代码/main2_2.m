% 最终推荐方案的确定程序二
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
SB = 100;
UB = 220;
ZB = UB.^2 ./ SB;

% 清理无用变量
clear S1_min S2_min S3_min S1_max S2_max S3_max SA_max

% ===================== 读取输电线参数 =====================
SDX = readmatrix("输电线参数.txt");
R_Plan5 = SDX(1:5,1:5);
X_Plan5 = SDX(1:5,6:10);
B_Plan5 = SDX(1:5,11:15);
R_Plan10 = SDX(6:10,1:5);
X_Plan10 = SDX(6:10,6:10);
B_Plan10 = SDX(6:10,11:15);
clear SDX

% ===================== 方案 5 潮流计算 =====================
disp('方案 5 潮流:')
SA_case_list = [S_A; 0.5*S_A];
qk = 1;
SA_case = SA_case_list(qk);
% A发电机出力（P扣厂用电，Q不扣）
PdA = real(SA_case) * (1 - k3_A);
QdA = imag(SA_case);

[U5,theta5,P5,Q5,Branch_Loss5] = CLJS_Plan5( ...
    Pd1, Pd2, Pd3,  Qd1, Qd2, Qd3, ...
    ZB, UB, SB, ...
    PdA, QdA, PdA_load, QdA_load, ...
    R_Plan5,X_Plan5,B_Plan5);
cosphi5 = P5./abs(P5+1j*Q5);
TABLE5 = [U5*UB,theta5,P5,Q5,cosphi5];
disp(TABLE5)

% ===================== 方案 10 潮流计算 =====================
disp('方案 10 潮流:')
SA_case_list = [S_A; 0.5*S_A];
qk = 1;
SA_case = SA_case_list(qk);
% A发电机出力（P扣厂用电，Q不扣）
PdA = real(SA_case) * (1 - k3_A);
QdA = imag(SA_case);
[U10,theta10,P10,Q10,Branch_Loss10] = CLJS_Plan10( ...
    Pd1, Pd2, Pd3,  Qd1, Qd2, Qd3, ...
    ZB, UB, SB, ...
    PdA, QdA, PdA_load, QdA_load, ...
    R_Plan10,X_Plan10,B_Plan10);
cosphi10 = P10./abs(P10+1j*Q10);
TABLE10 = [U10*UB,theta10,P10,Q10,cosphi10];
disp(TABLE10)

% ===================== 损耗与经济指标计算 =====================
% 支路有功损耗
dP5 = Branch_Loss5(:,3);
dP10 = Branch_Loss10(:,3);

% 年利用小时数
tao_5 = [4100 4100 3400 3400 3500 3500];
tao_10 = [4100 3400 4100 3400 3400 3400 3500];
dW5 = tao_5*dP5;
dW10 = tao_10*dP10;

% 资金回收系数计算
ii = 0.1;
nn = 25;
changshu = (ii*(1+ii)^nn)/((1+ii)^nn-1);

% -------------------- 方案 5 经济与电压指标 --------------------
f = zeros(number_of_point);
a = ones(number_of_point);
f(1,4) = 1; f(1,5) = 1;
f(2,4) = 1; f(2,5) = 1;
f(3,4) = 1; f(3,5) = 1;
a(1,4) = 2;
a(2,4) = 2;
f = f + f.';
a = a + a.';
zeta = 0.2;
kxi = zeros(number_of_point);
kxi(1,4) = 60; kxi(1,5) = 60;
kxi(2,4) = 60; kxi(3,4) = 60; kxi(2,5) = 60;
kxi = 10*(kxi+kxi.');
Im5 = sum(d.*f.*kxi.*(1+(a-1)*zeta),"all");
alpha = 0.5;
alpha1 = 0.05;
alpha2 = 0.04;
Cm5 = (alpha1+alpha2)*Im5 + alpha*dW5;
ACm5 = Cm5 + Im5*changshu;
dUmax5 = 100*max(abs((U5*UB-220))/220);

% -------------------- 方案 10 经济与电压指标 --------------------
f = zeros(number_of_point);
a = ones(number_of_point);
f(1,2) = 1; f(1,4) = 1; f(1,5) = 1;
f(2,4) = 1; f(2,5) = 1;
f(2,3) = 1;
f(3,4) = 1;
a(1,4) = 1; a(1,5) = 2;
a(2,4) = 1; a(2,5) = 2;
f = f + f.';
a = a + a.';
zeta = 0.2;
kxi = zeros(number_of_point);
kxi(1,2) = 60; kxi(1,4) = 60; kxi(1,5) = 60;
kxi(2,4) = 60; kxi(2,5) = 60; kxi(2,3) = 60;
kxi(3,4) = 60;
kxi = (kxi+kxi.')*10;
Im10 = sum(d.*f.*kxi.*(1+(a-1)*zeta),"all");
alpha = 0.5;
alpha1 = 0.05;
alpha2 = 0.04;
Cm10 = (alpha1+alpha2)*Im10 + alpha*dW10;
ACm10 = Cm10 + Im10*changshu;
dUmax10 = 100*max(abs((U10*UB-220))/220);

% ===================== 技术经济指标绘图 =====================
% figure
% yyaxis left
% bar([2,5,7,8],[dUmax5 dUmax10 0 0])
% title('方案技术与经济指标')
% ylabel('最大电压偏差量')
% ylim([0,13])
% xticks(1:7)
% labels = {'','','方案五','','','方案十'};
% xticklabels(labels)
% grid on;
% set(gca,'LineWidth',1.65,'FontSize',13.5);
% 
% yyaxis right
% bar([8,9,3,6,],[0 0 ACm5 ACm10])
% xlim([1,7])
% ylabel('年费用(千元)')
% legend('最大电压偏差量','年费用(千元)')

figure

yyaxis left
bar([2,5,7,8],[dUmax5 dUmax10 0 0])
title('方案技术与经济指标')
ylabel('最大电压偏差量')
ylim('auto')          % 左轴自适应
xticks([2.5 5.5])                 % 两组柱子的中心
xticklabels({'方案五','方案十'})
grid on;
set(gca,'LineWidth',1.65,'FontSize',13.5);

yyaxis right
bar([8,9,3,6],[0 0 ACm5 ACm10])
ylabel('年费用(千元)')
ylim('auto')          % 右轴自适应

xlim([1,7])
legend('最大电压偏差量','年费用(千元)')



figure

yyaxis left
bar(3, dUmax5)
title('方案五技术与经济指标')
ylabel('最大电压偏差量')
ylim('auto')

xticks(3)
xticklabels({'方案五'})
grid on;
set(gca,'LineWidth',1.65,'FontSize',13.5);

yyaxis right
bar(4, ACm5)
ylabel('年费用(千元)')
ylim('auto')

xlim([2,5])
legend('最大电压偏差量','年费用(千元)')




figure

yyaxis left
bar(3, dUmax10)
title('方案十技术与经济指标')
ylabel('最大电压偏差量')
ylim('auto')

xticks(3)
xticklabels({'方案十'})
grid on;
set(gca,'LineWidth',1.65,'FontSize',13.5);

yyaxis right
bar(4, ACm10)
ylabel('年费用(千元)')
ylim('auto')

xlim([2,5])
legend('最大电压偏差量','年费用(千元)')

% ===================== PCA 综合评价 =====================
% 指标归一化
gjU = [dUmax5;dUmax10];
gjAC = [ACm5;ACm10];
gjU = 1-(gjU-min(gjU))/(max(gjU)-min(gjU));
gjAC = 1-(gjAC-min(gjAC))/(max(gjAC)-min(gjAC));

% 主成分分析
gj = [gjU,gjAC];
gj = zscore(gj);
r = corrcoef(gj);
[x,y,z] = pcacov(r);
F = repmat(sign(sum(x)),size(x,1),1);
x = x.*F;
num = 2;
df = gj*x(:,[1,num]);
tf = df*z(1:num)/100'; % 综合得分
[stf,ind] = sort(tf,'descend'); % 从高到低排序
stf
ind

% PCA 结果绘图
figure
bar(1,stf(1),'g')
hold on
bar(2,stf(2),'m')
xticks(1:2)
xticklabels({'方案五','方案十'})
ylabel('综合得分')
grid on;
set(gca,'LineWidth',1.65,'FontSize',13.5);

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数：潮流计算核心函数
% 输入：mpc - 潮流计算结构体
% 输出：U - 节点电压幅值，theta - 节点电压相角
%       P - 节点注入有功功率，Q - 节点注入无功功率
%       Branch_Loss - 支路损耗矩阵（起点/终点/有功损/无功损/输送有功/输送无功）
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
% 函数：方案5潮流计算封装（含导线参数修正）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U,theta,P,Q,Branch_Loss] = CLJS_Plan5( ...
    Pd1, Pd2, Pd3,  Qd1, Qd2, Qd3, ...
    ZB, UB, SB, ...
    PdA, QdA, PdA_load, QdA_load, ...
    R_Plan5,X_Plan5,B_Plan5)
% 阻抗导纳归算到基准值
r = R_Plan5./ZB;
x = X_Plan5./ZB;
b = B_Plan5.*ZB; % YB=1/ZB
% 修正2回线路参数
r(1,4) = r(1,4)/2; x(1,4) = x(1,4)/2; b(1,4) = b(1,4)*2;
r(2,4) = r(2,4)/2; x(2,4) = x(2,4)/2; b(2,4) = b(2,4)*2;

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
    1 4 r(1,4) x(1,4) b(1,4) 0 0 0 0 0 1 -360 360;  % 支路1-2（匹配f(1,2)=1）
    1 5 r(1,5) x(1,5) b(1,5) 0 0 0 0 0 1 -360 360;
    2 4 r(2,4) x(2,4) b(2,4) 0 0 0 0 0 1 -360 360;  % 支路2-4（匹配f(2,4)=1）
    2 5 r(2,5) x(2,5) b(2,5) 0 0 0 0 0 1 -360 360;  % 支路2-5（匹配f(2,5)=1）
    3 4 r(3,4) x(3,4) b(3,4) 0 0 0 0 0 1 -360 360;  % 支路3-4（匹配f(3,4)=1）
    3 5 r(3,5) x(3,5) b(3,5) 0 0 0 0 0 1 -360 360;
    ];

% 调用潮流计算函数
[U,theta,P,Q,Branch_Loss] = Power_Flow_Calculation(mpc);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数：方案10潮流计算封装（含导线参数修正）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U,theta,P,Q,Branch_Loss] = CLJS_Plan10( ...
    Pd1, Pd2, Pd3,  Qd1, Qd2, Qd3, ...
    ZB, UB, SB, ...
    PdA, QdA, PdA_load, QdA_load, ...
    R_Plan10,X_Plan10,B_Plan10)
    % 阻抗导纳归算到基准值
    r = R_Plan10./ZB;
    x = X_Plan10./ZB;
    b = B_Plan10.*ZB; % YB=1/ZB
    % 修正2回线路参数
    % r(1,4) = r(1,4)/2; x(1,4) = x(1,4)/2; b(1,4) = b(1,4)*2;
    r(1,5) = r(1,5)/2; x(1,5) = x(1,5)/2; b(1,5) = b(1,5)*2;
    % r(2,4) = r(2,4)/2; x(2,4) = x(2,4)/2; b(2,4) = b(2,4)*2;
    r(2,5) = r(2,5)/2; x(2,5) = x(2,5)/2; b(2,5) = b(2,5)*2;
    
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
    1 2 r(1,2) x(1,2) b(1,2) 0 0 0 0 0 1 -360 360;  % 支路1-2（匹配f(1,2)=1）
    1 4 r(1,4) x(1,4) b(1,4) 0 0 0 0 0 1 -360 360;  
    1 5 r(1,5) x(1,5) b(1,5) 0 0 0 0 0 1 -360 360;  
    2 3 r(2,3) x(2,3) b(2,3) 0 0 0 0 0 1 -360 360;
    2 4 r(2,4) x(2,4) b(2,4) 0 0 0 0 0 1 -360 360;  % 支路2-4（匹配f(2,4)=1）
    2 5 r(2,5) x(2,5) b(2,5) 0 0 0 0 0 1 -360 360;  % 支路2-5（匹配f(2,5)=1）
    3 4 r(3,4) x(3,4) b(3,4) 0 0 0 0 0 1 -360 360;
    ];

    
    % 调用潮流计算函数
    [U,theta,P,Q,Branch_Loss] = Power_Flow_Calculation(mpc);
end