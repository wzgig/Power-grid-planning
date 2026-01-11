% 初步接线方案的选择程序
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

%% =========================
%  结果表格初始化（7个方案，3类指标）
%% =========================
Table = zeros(10,3);

% 仅两种工况权重（A两台满发 / A一台满发）
alpha = [0.7; 0.3];

%%% 方案 1 计算
DNZL = zeros(2,1);
XLSH = zeros(2,1);
XLTZ = zeros(2,1);
% 两种工况下 A 的"视在功率总额"（两台 / 一台）
SA_case_list = [S_A; 0.5*S_A];

for qk = 1:2
    SA_case = SA_case_list(qk);

    % A发电机出力（P扣厂用电，Q不扣）
    PdA = real(SA_case) * (1 - k3_A);  % 这里沿用你 Plan_1 参数名 PdA/QdA
    QdA = imag(SA_case);

    % 调用 Plan_1（Slack负荷置0）
    [dnzl, xltz, xlsh] = Plan_1( ...
        Pd1, Pd2, PdS,  Qd1, Qd2, QdS, ...
        Z_all, ZB, UB, SB, ...
        PdA, QdA, PdA_load, QdA_load, ...
        number_of_point, d);

    DNZL(qk) = dnzl;
    XLTZ(qk) = xltz;
    XLSH(qk) = xlsh;
end

Table(1,1) = dot(alpha,DNZL);
Table(1,2) = dot(alpha,XLTZ);
Table(1,3) = dot(alpha,XLSH);


%%% 方案 2 计算
DNZL = zeros(2,1);
XLSH = zeros(2,1);
XLTZ = zeros(2,1);
% 两种工况下 A 的"视在功率总额"（两台 / 一台）
SA_case_list = [S_A; 0.5*S_A];

for qk = 1:2
    SA_case = SA_case_list(qk);

    % A发电机出力（P扣厂用电，Q不扣）
    PdA = real(SA_case) * (1 - k3_A);  % 这里沿用你 Plan_1 参数名 PdA/QdA
    QdA = imag(SA_case);

    % 调用 Plan_1（Slack负荷置0）
    [dnzl, xltz, xlsh] = Plan_2( ...
        Pd1, Pd2, PdS,  Qd1, Qd2, QdS, ...
        Z_all, ZB, UB, SB, ...
        PdA, QdA, PdA_load, QdA_load, ...
        number_of_point, d);

    DNZL(qk) = dnzl;
    XLTZ(qk) = xltz;
    XLSH(qk) = xlsh;
end

Table(2,1) = dot(alpha,DNZL);
Table(2,2) = dot(alpha,XLTZ);
Table(2,3) = dot(alpha,XLSH);

%%% 方案 3 计算
DNZL = zeros(2,1);
XLSH = zeros(2,1);
XLTZ = zeros(2,1);
% 两种工况下 A 的"视在功率总额"（两台 / 一台）
SA_case_list = [S_A; 0.5*S_A];

for qk = 1:2
    SA_case = SA_case_list(qk);

    % A发电机出力（P扣厂用电，Q不扣）
    PdA = real(SA_case) * (1 - k3_A);  % 这里沿用你 Plan_1 参数名 PdA/QdA
    QdA = imag(SA_case);

    % 调用 Plan_1（Slack负荷置0）
    [dnzl, xltz, xlsh] = Plan_3( ...
        Pd1, Pd2, PdS,  Qd1, Qd2, QdS, ...
        Z_all, ZB, UB, SB, ...
        PdA, QdA, PdA_load, QdA_load, ...
        number_of_point, d);

    DNZL(qk) = dnzl;
    XLTZ(qk) = xltz;
    XLSH(qk) = xlsh;
end

Table(3,1) = dot(alpha,DNZL);
Table(3,2) = dot(alpha,XLTZ);
Table(3,3) = dot(alpha,XLSH);


%%% 方案 4 计算
DNZL = zeros(2,1);
XLSH = zeros(2,1);
XLTZ = zeros(2,1);
% 两种工况下 A 的"视在功率总额"（两台 / 一台）
SA_case_list = [S_A; 0.5*S_A];

for qk = 1:2
    SA_case = SA_case_list(qk);

    % A发电机出力（P扣厂用电，Q不扣）
    PdA = real(SA_case) * (1 - k3_A);  % 这里沿用你 Plan_1 参数名 PdA/QdA
    QdA = imag(SA_case);

    % 调用 Plan_1（Slack负荷置0）
    [dnzl, xltz, xlsh] = Plan_4( ...
        Pd1, Pd2, PdS,  Qd1, Qd2, QdS, ...
        Z_all, ZB, UB, SB, ...
        PdA, QdA, PdA_load, QdA_load, ...
        number_of_point, d);

    DNZL(qk) = dnzl;
    XLTZ(qk) = xltz;
    XLSH(qk) = xlsh;
end
Table(4,1) = dot(alpha,DNZL);
Table(4,2) = dot(alpha,XLTZ);
Table(4,3) = dot(alpha,XLSH);

%%% 方案 5 计算
DNZL = zeros(2,1);
XLSH = zeros(2,1);
XLTZ = zeros(2,1);
% 两种工况下 A 的"视在功率总额"（两台 / 一台）
SA_case_list = [S_A; 0.5*S_A];

for qk = 1:2
    SA_case = SA_case_list(qk);

    % A发电机出力（P扣厂用电，Q不扣）
    PdA = real(SA_case) * (1 - k3_A);  % 这里沿用你 Plan_1 参数名 PdA/QdA
    QdA = imag(SA_case);

    % 调用 Plan_1（Slack负荷置0）
    [dnzl, xltz, xlsh] = Plan_5( ...
        Pd1, Pd2, PdS,  Qd1, Qd2, QdS, ...
        Z_all, ZB, UB, SB, ...
        PdA, QdA, PdA_load, QdA_load, ...
        number_of_point, d);

    DNZL(qk) = dnzl;
    XLTZ(qk) = xltz;
    XLSH(qk) = xlsh;
end
Table(5,1) = dot(alpha,DNZL);
Table(5,2) = dot(alpha,XLTZ);
Table(5,3) = dot(alpha,XLSH);

%%% 方案 6 计算
DNZL = zeros(2,1);
XLSH = zeros(2,1);
XLTZ = zeros(2,1);
% 两种工况下 A 的"视在功率总额"（两台 / 一台）
SA_case_list = [S_A; 0.5*S_A];

for qk = 1:2
    SA_case = SA_case_list(qk);

    % A发电机出力（P扣厂用电，Q不扣）
    PdA = real(SA_case) * (1 - k3_A);  % 这里沿用你 Plan_1 参数名 PdA/QdA
    QdA = imag(SA_case);

    % 调用 Plan_1（Slack负荷置0）
    [dnzl, xltz, xlsh] = Plan_6( ...
        Pd1, Pd2, PdS,  Qd1, Qd2, QdS, ...
        Z_all, ZB, UB, SB, ...
        PdA, QdA, PdA_load, QdA_load, ...
        number_of_point, d);

    DNZL(qk) = dnzl;
    XLTZ(qk) = xltz;
    XLSH(qk) = xlsh;
end
Table(6,1) = dot(alpha,DNZL);
Table(6,2) = dot(alpha,XLTZ);
Table(6,3) = dot(alpha,XLSH);

%%% 方案 7 计算
DNZL = zeros(2,1);
XLSH = zeros(2,1);
XLTZ = zeros(2,1);
% 两种工况下 A 的"视在功率总额"（两台 / 一台）
SA_case_list = [S_A; 0.5*S_A];

for qk = 1:2
    SA_case = SA_case_list(qk);

    % A发电机出力（P扣厂用电，Q不扣）
    PdA = real(SA_case) * (1 - k3_A);  % 这里沿用你 Plan_1 参数名 PdA/QdA
    QdA = imag(SA_case);

    % 调用 Plan_1（Slack负荷置0）
    [dnzl, xltz, xlsh] = Plan_7( ...
        Pd1, Pd2, PdS,  Qd1, Qd2, QdS, ...
        Z_all, ZB, UB, SB, ...
        PdA, QdA, PdA_load, QdA_load, ...
        number_of_point, d);

    DNZL(qk) = dnzl;
    XLTZ(qk) = xltz;
    XLSH(qk) = xlsh;
end
Table(7,1) = dot(alpha,DNZL);
Table(7,2) = dot(alpha,XLTZ);
Table(7,3) = dot(alpha,XLSH);

%%% 方案 8 计算
DNZL = zeros(2,1);
XLSH = zeros(2,1);
XLTZ = zeros(2,1);
% 两种工况下 A 的"视在功率总额"（两台 / 一台）
SA_case_list = [S_A; 0.5*S_A];

for qk = 1:2
    SA_case = SA_case_list(qk);

    % A发电机出力（P扣厂用电，Q不扣）
    PdA = real(SA_case) * (1 - k3_A);  % 这里沿用你 Plan_1 参数名 PdA/QdA
    QdA = imag(SA_case);

    % 调用 Plan_1（Slack负荷置0）
    [dnzl, xltz, xlsh] = Plan_8( ...
        Pd1, Pd2, PdS,  Qd1, Qd2, QdS, ...
        Z_all, ZB, UB, SB, ...
        PdA, QdA, PdA_load, QdA_load, ...
        number_of_point, d);

    DNZL(qk) = dnzl;
    XLTZ(qk) = xltz;
    XLSH(qk) = xlsh;
end
Table(8,1) = dot(alpha,DNZL);
Table(8,2) = dot(alpha,XLTZ);
Table(8,3) = dot(alpha,XLSH);

%%% 方案 9 计算
DNZL = zeros(2,1);
XLSH = zeros(2,1);
XLTZ = zeros(2,1);
% 两种工况下 A 的"视在功率总额"（两台 / 一台）
SA_case_list = [S_A; 0.5*S_A];

for qk = 1:2
    SA_case = SA_case_list(qk);

    % A发电机出力（P扣厂用电，Q不扣）
    PdA = real(SA_case) * (1 - k3_A);  % 这里沿用你 Plan_1 参数名 PdA/QdA
    QdA = imag(SA_case);

    % 调用 Plan_1（Slack负荷置0）
    [dnzl, xltz, xlsh] = Plan_9( ...
        Pd1, Pd2, PdS,  Qd1, Qd2, QdS, ...
        Z_all, ZB, UB, SB, ...
        PdA, QdA, PdA_load, QdA_load, ...
        number_of_point, d);

    DNZL(qk) = dnzl;
    XLTZ(qk) = xltz;
    XLSH(qk) = xlsh;
end
Table(9,1) = dot(alpha,DNZL);
Table(9,2) = dot(alpha,XLTZ);
Table(9,3) = dot(alpha,XLSH);


%%% 方案 10 计算
DNZL = zeros(2,1);
XLSH = zeros(2,1);
XLTZ = zeros(2,1);
% 两种工况下 A 的"视在功率总额"（两台 / 一台）
SA_case_list = [S_A; 0.5*S_A];

for qk = 1:2
    SA_case = SA_case_list(qk);

    % A发电机出力（P扣厂用电，Q不扣）
    PdA = real(SA_case) * (1 - k3_A);  % 这里沿用你 Plan_1 参数名 PdA/QdA
    QdA = imag(SA_case);

    % 调用 Plan_1（Slack负荷置0）
    [dnzl, xltz, xlsh] = Plan_10( ...
        Pd1, Pd2, PdS,  Qd1, Qd2, QdS, ...
        Z_all, ZB, UB, SB, ...
        PdA, QdA, PdA_load, QdA_load, ...
        number_of_point, d);

    DNZL(qk) = dnzl;
    XLTZ(qk) = xltz;
    XLSH(qk) = xlsh;
end
Table(10,1) = dot(alpha,DNZL);
Table(10,2) = dot(alpha,XLTZ);
Table(10,3) = dot(alpha,XLSH);

% 跳过方案4
figure
x = (1:10).';
x_filter = x([1:3,5:10]);       % 过滤后的方案编号
table_filter = Table([1:3,5:10],1);  % 过滤后的技术指标数据
plot(x_filter, table_filter, '-ob', 'LineWidth', 1.5)  % 仅画1-3、5-10号方案
xticks(1:10)  % 刻度仍显示1-10（保持坐标轴完整）
labels = {'1','2','3','4','5','6','7','8','9','10'};
xticklabels(labels)
grid on
box on
title('技术指标')
set(gca,'LineWidth',1.15,'FontSize',13);
xlabel('方案编号')
ylabel('电能质量(10^9·mVA)')

%% 先定义公共参数（避免重复）
x = (1:10).';  % 原始10个方案编号
labels = {'1','2','3','4','5','6','7','8','9','10'};  % 完整刻度标签
% 核心：定义过滤索引（跳过第四个点，保留1-3、5-10）
filter_idx = [1,2,3,5:10];  
x_filter = x(filter_idx);  % 过滤后的方案编号（无方案4）

%% 经济指标图（双轴：线路投资+线路损耗，跳过方案4）
figure
% ===== 左侧y轴：线路投资（跳过方案4） =====
yyaxis left
% 绘制1-3、5-10号方案的线路投资（剔除方案4）
plot(x_filter, Table(filter_idx,2), '-o', "LineWidth", 1.5)  
title('经济指标')
xlabel('方案编号')
ylabel('线路投资(km)')
xticks(1:10)                                 % 刻度仍显示1-10（坐标轴完整）
xticklabels(labels)
grid on;
set(gca, 'LineWidth', 1.5, 'FontSize', 13.5);  

% ===== 右侧y轴：线路损耗（跳过方案4） =====
yyaxis right
% 绘制1-3、5-10号方案的线路损耗（剔除方案4）
plot(x_filter, Table(filter_idx,3), '-', "LineWidth", 1.5, "Marker", "diamond")  
ylabel('线路损耗(10^{15}·mV^2A^2)')  
set(gca, 'LineWidth', 1.5, 'FontSize', 13.5);  

%% 第二个经济指标图（双轴：跳过方案4，保持逻辑一致）
figure
% ===== 左侧y轴：线路投资（跳过方案4） =====
yyaxis left
plot(x_filter, Table(filter_idx,2), '-o', "LineWidth", 1.5)  
xlabel('方案编号')
xticks(1:10)                                
xticklabels(labels)                         
grid on;
set(gca, 'LineWidth', 1.5, 'FontSize', 13.5);  

% ===== 右侧y轴：线路损耗 =====
yyaxis right
plot(x_filter, Table(filter_idx,3), '-', "LineWidth", 1.5, "Marker", "diamond")  
ylabel('线路损耗(10^{15}·mV^2A^2)')
set(gca, 'LineWidth', 1.5, 'FontSize', 13.5);

% 候选方案：跳过4
candidate_idx = setdiff(1:size(Table,1), 4);

% 只取候选数据
T = Table(candidate_idx, :);     % T(:,1)=power, T(:,2)=investment, T(:,3)=loss

% -------- 技术指标：最小（只在候选里比）
[DNZL_MIN, loc] = min(T(:,1));
IND = candidate_idx(loc);        % 映射回原表行号

% -------- 经济指标：归一化（只在候选里算min/max）
power = (T(:,1) - min(T(:,1))) / (max(T(:,1)) - min(T(:,1)))
investment   = (T(:,2) - min(T(:,2))) / (max(T(:,2)) - min(T(:,2)))
loss  = (T(:,3) - min(T(:,3))) / (max(T(:,3)) - min(T(:,3)))

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 辅助函数：参数计算函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dnzl,xltz,xlsh] = canshujisuan(d,f,a,Branch_Loss,number_of_point)
S = zeros(number_of_point);
qidian = Branch_Loss(:,1);
zhongdian = Branch_Loss(:,2);
n = length(qidian);
for i = 1:n
    S(qidian(i),zhongdian(i)) = ...
        sqrt((Branch_Loss(i,3)).^2+(Branch_Loss(i,4)).^2);
end
dnzl = max(max(d.*f.*S./a));
xlsh = sum(d.*f.*S.*S./a,'all');
xltz = sum(d.*f.*a,'all')/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 辅助函数：潮流计算函数
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

% 支路损耗： 起点,终点,有功网损， 无功网损
Branch_Loss = [
    Power_flow_result.branch(:,1),Power_flow_result.branch(:,2),...
    Power_flow_result.branch(:,14)+Power_flow_result.branch(:,16),...
    Power_flow_result.branch(:,15)+Power_flow_result.branch(:,17)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 方案1计算函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dnzl, xltz, xlsh] = Plan_1( ...
    Pd1, Pd2, Pd3,  ...
    Qd1, Qd2, Qd3,  ...
    Z_all, ZB, UB, SB, ...
    PdA, QdA, ...
    PdA_load, QdA_load, ...
    number_of_point, d)

f = zeros(number_of_point); % 节点数×节点数的全0矩阵
a = ones(number_of_point);% 全1矩阵，节点之间是否双回路接线
% 定义支路连接：f(i,j)=1表示节点i-j有支路
f(1,2) = 1; f(1,4) = 1; f(1,5) = 1;
f(2,4) = 1; f(2,5) = 1;
f(3,4) = 1;
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
    1 4 r(1,4) x(1,4) 0 0 0 0 0 0 1 -360 360;  % 支路1-4（匹配f(1,4)=1）
    1 5 r(1,5) x(1,5) 0 0 0 0 0 0 1 -360 360;  % 支路1-5（匹配f(1,5)=1）
    2 4 r(2,4) x(2,4) 0 0 0 0 0 0 1 -360 360;  % 支路2-4（匹配f(2,4)=1）
    2 5 r(2,5) x(2,5) 0 0 0 0 0 0 1 -360 360;  % 支路2-5（匹配f(2,5)=1）
    3 4 r(3,4) x(3,4) 0 0 0 0 0 0 1 -360 360;  % 支路3-4（匹配f(3,4)=1）
    ];

% 潮流计算与参数计算
[~,~,~,~,Branch_Loss] = Power_Flow_Calculation(mpc);
[dnzl,xltz,xlsh] = canshujisuan(d,f,a,Branch_Loss,number_of_point);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 方案2计算函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dnzl, xltz, xlsh] = Plan_2( ...
    Pd1, Pd2, Pd3,  ...
    Qd1, Qd2, Qd3,  ...
    Z_all, ZB, UB, SB, ...
    PdA, QdA, ...
    PdA_load, QdA_load, ...
    number_of_point, d)

f = zeros(number_of_point); % 节点数×节点数的全0矩阵
a = ones(number_of_point);% 全1矩阵，节点之间是否双回路接线
% 定义支路连接：f(i,j)=1表示节点i-j有支路
f(1,2) = 1; f(1,4) = 1;
f(2,4) = 1; f(2,5) = 1;
f(3,4) = 1;
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
    1 4 r(1,4) x(1,4) 0 0 0 0 0 0 1 -360 360;  % 支路1-4（匹配f(1,4)=1）
    2 4 r(2,4) x(2,4) 0 0 0 0 0 0 1 -360 360;  % 支路2-4（匹配f(2,4)=1）
    2 5 r(2,5) x(2,5) 0 0 0 0 0 0 1 -360 360;  % 支路2-5（匹配f(2,5)=1）
    3 4 r(3,4) x(3,4) 0 0 0 0 0 0 1 -360 360;  % 支路3-4（匹配f(3,4)=1）
    ];

% 潮流计算与参数计算
[~,~,~,~,Branch_Loss] = Power_Flow_Calculation(mpc);
[dnzl,xltz,xlsh] = canshujisuan(d,f,a,Branch_Loss,number_of_point);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 方案3计算函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dnzl, xltz, xlsh] = Plan_3( ...
    Pd1, Pd2, Pd3,  ...
    Qd1, Qd2, Qd3,  ...
    Z_all, ZB, UB, SB, ...
    PdA, QdA, ...
    PdA_load, QdA_load, ...
    number_of_point, d)

f = zeros(number_of_point); % 节点数×节点数的全0矩阵
a = ones(number_of_point);% 全1矩阵，节点之间是否双回路接线
% 定义支路连接：f(i,j)=1表示节点i-j有支路
f(1,4) = 1; f(1,5) = 1;
f(2,4) = 1; f(2,5) = 1;
f(3,4) = 1;
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
    1 5 r(1,5) x(1,5) 0 0 0 0 0 0 1 -360 360;  % 支路1-4（匹配f(1,4)=1）
    2 4 r(2,4) x(2,4) 0 0 0 0 0 0 1 -360 360;  % 支路2-4（匹配f(2,4)=1）
    2 5 r(2,5) x(2,5) 0 0 0 0 0 0 1 -360 360;  % 支路2-5（匹配f(2,5)=1）
    3 4 r(3,4) x(3,4) 0 0 0 0 0 0 1 -360 360;  % 支路3-4（匹配f(3,4)=1）
    ];

% 潮流计算与参数计算
[~,~,~,~,Branch_Loss] = Power_Flow_Calculation(mpc);
[dnzl,xltz,xlsh] = canshujisuan(d,f,a,Branch_Loss,number_of_point);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 方案4计算函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dnzl, xltz, xlsh] = Plan_4( ...
    Pd1, Pd2, Pd3,  ...
    Qd1, Qd2, Qd3,  ...
    Z_all, ZB, UB, SB, ...
    PdA, QdA, ...
    PdA_load, QdA_load, ...
    number_of_point, d)

f = zeros(number_of_point); % 节点数×节点数的全0矩阵
a = ones(number_of_point);% 全1矩阵，节点之间是否双回路接线
% 定义支路连接：f(i,j)=1表示节点i-j有支路
f(1,4) = 1;
f(2,4) = 1; f(2,5) = 1;
f(3,4) = 1;
a(1,4) = 2;
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
    2 4 r(2,4) x(2,4) 0 0 0 0 0 0 1 -360 360;  % 支路2-4（匹配f(2,4)=1）
    2 5 r(2,5) x(2,5) 0 0 0 0 0 0 1 -360 360;  % 支路2-5（匹配f(2,5)=1）
    3 4 r(3,4) x(3,4) 0 0 0 0 0 0 1 -360 360;  % 支路3-4（匹配f(3,4)=1）
    ];

% 潮流计算与参数计算
[~,~,~,~,Branch_Loss] = Power_Flow_Calculation(mpc);
[dnzl,xltz,xlsh] = canshujisuan(d,f,a,Branch_Loss,number_of_point);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 方案5计算函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dnzl, xltz, xlsh] = Plan_5( ...
    Pd1, Pd2, Pd3,  ...
    Qd1, Qd2, Qd3,  ...
    Z_all, ZB, UB, SB, ...
    PdA, QdA, ...
    PdA_load, QdA_load, ...
    number_of_point, d)

f = zeros(number_of_point); % 节点数×节点数的全0矩阵
a = ones(number_of_point);% 全1矩阵，节点之间是否双回路接线
% 定义支路连接：f(i,j)=1表示节点i-j有支路
f(1,4) = 1; f(1,5) = 1;
f(2,4) = 1; f(2,5) = 1;
f(3,4) = 1;
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
    ];

% 潮流计算与参数计算
[~,~,~,~,Branch_Loss] = Power_Flow_Calculation(mpc);
[dnzl,xltz,xlsh] = canshujisuan(d,f,a,Branch_Loss,number_of_point);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 方案6计算函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dnzl, xltz, xlsh] = Plan_6( ...
    Pd1, Pd2, Pd3,  ...
    Qd1, Qd2, Qd3,  ...
    Z_all, ZB, UB, SB, ...
    PdA, QdA, ...
    PdA_load, QdA_load, ...
    number_of_point, d)

f = zeros(number_of_point); % 节点数×节点数的全0矩阵
a = ones(number_of_point);% 全1矩阵，节点之间是否双回路接线
% 定义支路连接：f(i,j)=1表示节点i-j有支路
f(1,2) = 1; f(1,4) = 1;
f(2,4) = 1; f(2,5) = 1;
f(3,4) = 1;
a(1,4) = 2;
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
    2 4 r(2,4) x(2,4) 0 0 0 0 0 0 1 -360 360;  % 支路2-4（匹配f(2,4)=1）
    2 5 r(2,5) x(2,5) 0 0 0 0 0 0 1 -360 360;  % 支路2-5（匹配f(2,5)=1）
    3 4 r(3,4) x(3,4) 0 0 0 0 0 0 1 -360 360;  % 支路3-4（匹配f(3,4)=1）
    ];

% 潮流计算与参数计算
[~,~,~,~,Branch_Loss] = Power_Flow_Calculation(mpc);
[dnzl,xltz,xlsh] = canshujisuan(d,f,a,Branch_Loss,number_of_point);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 方案7计算函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dnzl, xltz, xlsh] = Plan_7( ...
    Pd1, Pd2, Pd3,  ...
    Qd1, Qd2, Qd3,  ...
    Z_all, ZB, UB, SB, ...
    PdA, QdA, ...
    PdA_load, QdA_load, ...
    number_of_point, d)

f = zeros(number_of_point); % 节点数×节点数的全0矩阵
a = ones(number_of_point);% 全1矩阵，节点之间是否双回路接线
% 定义支路连接：f(i,j)=1表示节点i-j有支路
f(1,2) = 1; f(1,4) = 1; f(1,5) = 1;
f(2,4) = 1; f(2,5) = 1;
f(3,4) = 1;
a(1,2) = 2;
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
    2 4 r(2,4) x(2,4) 0 0 0 0 0 0 1 -360 360;  % 支路2-4（匹配f(2,4)=1）
    2 5 r(2,5) x(2,5) 0 0 0 0 0 0 1 -360 360;  % 支路2-5（匹配f(2,5)=1）
    3 4 r(3,4) x(3,4) 0 0 0 0 0 0 1 -360 360;  % 支路3-4（匹配f(3,4)=1）
    ];

% 潮流计算与参数计算
[~,~,~,~,Branch_Loss] = Power_Flow_Calculation(mpc);
[dnzl,xltz,xlsh] = canshujisuan(d,f,a,Branch_Loss,number_of_point);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 方案8计算函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dnzl, xltz, xlsh] = Plan_8( ...
    Pd1, Pd2, Pd3,  ...
    Qd1, Qd2, Qd3,  ...
    Z_all, ZB, UB, SB, ...
    PdA, QdA, ...
    PdA_load, QdA_load, ...
    number_of_point, d)

f = zeros(number_of_point); % 节点数×节点数的全0矩阵
a = ones(number_of_point);% 全1矩阵，节点之间是否双回路接线
% 定义支路连接：f(i,j)=1表示节点i-j有支路
f(1,2) = 1; f(1,4) = 1; f(1,5) = 1;
f(2,4) = 1; f(2,5) = 1;
f(2,3) = 1;
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
    1 2 r(1,2) x(1,2) 0 0 0 0 0 0 1 -360 360;  % 支路1-2（匹配f(1,2)=1）
    1 4 r(1,4) x(1,4) 0 0 0 0 0 0 1 -360 360;  
    1 5 r(1,5) x(1,5) 0 0 0 0 0 0 1 -360 360;  
    2 3 r(2,3) x(2,3) 0 0 0 0 0 0 1 -360 360;
    2 4 r(2,4) x(2,4) 0 0 0 0 0 0 1 -360 360;  % 支路2-4（匹配f(2,4)=1）
    2 5 r(2,5) x(2,5) 0 0 0 0 0 0 1 -360 360;  % 支路2-5（匹配f(2,5)=1）
    ];

% 潮流计算与参数计算
[~,~,~,~,Branch_Loss] = Power_Flow_Calculation(mpc);
[dnzl,xltz,xlsh] = canshujisuan(d,f,a,Branch_Loss,number_of_point);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 方案9计算函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dnzl, xltz, xlsh] = Plan_9( ...
    Pd1, Pd2, Pd3,  ...
    Qd1, Qd2, Qd3,  ...
    Z_all, ZB, UB, SB, ...
    PdA, QdA, ...
    PdA_load, QdA_load, ...
    number_of_point, d)

f = zeros(number_of_point); % 节点数×节点数的全0矩阵
a = ones(number_of_point);% 全1矩阵，节点之间是否双回路接线
% 定义支路连接：f(i,j)=1表示节点i-j有支路
f(1,2) = 1; f(1,4) = 1; f(1,5) = 1;
f(2,4) = 1; f(2,5) = 1;
f(2,3) = 1;
a(1,4) = 2; a(1,5) = 2;
a(2,4) = 2; a(2,5) = 2;
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
    ];

% 潮流计算与参数计算
[~,~,~,~,Branch_Loss] = Power_Flow_Calculation(mpc);
[dnzl,xltz,xlsh] = canshujisuan(d,f,a,Branch_Loss,number_of_point);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 方案10计算函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dnzl, xltz, xlsh] = Plan_10( ...
    Pd1, Pd2, Pd3,  ...
    Qd1, Qd2, Qd3,  ...
    Z_all, ZB, UB, SB, ...
    PdA, QdA, ...
    PdA_load, QdA_load, ...
    number_of_point, d)

f = zeros(number_of_point); % 节点数×节点数的全0矩阵
a = ones(number_of_point);% 全1矩阵，节点之间是否双回路接线
% 定义支路连接：f(i,j)=1表示节点i-j有支路
f(1,2) = 1; f(1,4) = 1; f(1,5) = 1;
f(2,4) = 1; f(2,5) = 1;
f(2,3) = 1;
f(3,4) = 1;
a(1,4) = 2; a(1,5) = 2;
a(2,4) = 2; a(2,5) = 2;
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
[~,~,~,~,Branch_Loss] = Power_Flow_Calculation(mpc);
[dnzl,xltz,xlsh] = canshujisuan(d,f,a,Branch_Loss,number_of_point);
end