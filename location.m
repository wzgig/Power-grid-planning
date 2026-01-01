%% 地理位置关系图（终极兼容版：打印坐标+修复所有已知问题）
clear; clc; close all;

% ===================== 1. 坐标计算（核心修改：筛选A的有效坐标） =====================
S = [0, 0]; % S为原点
node2 = [140, 0]; % S到2的距离为140

% 节点1坐标（S=160，2=130）
x1 = (160^2 - 130^2 + 140^2) / (2*140); 
y1 = sqrt(160^2 - x1^2); 
node1 = [x1, y1];

% 节点3坐标（对称）
node3 = [x1, -y1];

% 节点A坐标（符号求解+筛选有效解）
syms x y
eq1 = (x - 140)^2 + y^2 == 140^2;  
eq2 = (x - x1)^2 + (y - y1)^2 == 160^2; 
sol = solve([eq1, eq2], [x, y], 'Real', true);

% 将解转为数值数组，筛选掉原点附近的无效解
x_sol = double(sol.x);
y_sol = double(sol.y);
% 保留距离原点大于10的解（排除(0,0)无效解）
valid_idx = sqrt(x_sol.^2 + y_sol.^2) > 10;
A = [x_sol(valid_idx), y_sol(valid_idx)];

% ===================== 新增功能：打印每个点的坐标（格式清晰） =====================
disp('=====================================');
disp('           各节点坐标（保留6位小数）');
disp('=====================================');
% 格式化打印，避免科学计数法，保留6位小数
disp(['S（原点）坐标：(', num2str(S(1), '%.6f'), ', ', num2str(S(2), '%.6f'), ')']);
disp(['节点1      坐标：(', num2str(node1(1), '%.6f'), ', ', num2str(node1(2), '%.6f'), ')']);
disp(['节点2      坐标：(', num2str(node2(1), '%.6f'), ', ', num2str(node2(2), '%.6f'), ')']);
disp(['节点3      坐标：(', num2str(node3(1), '%.6f'), ', ', num2str(node3(2), '%.6f'), ')']);
disp(['节点A      坐标：(', num2str(A(1), '%.6f'), ', ', num2str(A(2), '%.6f'), ')']);
disp('=====================================');

% ===================== 2. 绘图环境配置（全版本兼容） =====================
figure('Color', [0.98, 0.98, 0.98]); % 浅灰背景（护眼）
ax = gca;
hold on; grid on; axis equal;

% 坐标轴美化（核心修复：移除错误属性，全版本兼容）
ax.LineWidth = 1; % 坐标轴线条粗细
ax.GridAlpha = 0.3; % 网格透明度（弱化网格）
ax.GridColor = [0.7, 0.7, 0.7]; % 网格浅灰色
ax.FontName = '微软雅黑'; % 统一设置坐标轴所有文字（刻度+标签）的字体
ax.FontSize = 10; % 刻度标签字体大小
ax.Box = 'off'; % 隐藏坐标轴外框

% 轴标签（单独设置字号/加粗，继承FontName）
xlabel('X 坐标 (km)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y 坐标 (km)', 'FontSize', 12, 'FontWeight', 'bold');

% 标题美化
title('发电厂、变电所相对地理位置关系图', ...
    'FontName', '微软雅黑', 'FontSize', 14, 'FontWeight', 'bold', ...
    'Color', [0.1, 0.2, 0.4]); % 深蓝色标题

% ===================== 3. 绘制节点（立体+配色优化） =====================
% 配色方案（专业莫兰迪色，避免刺眼）
color_S = [0.2, 0.4, 0.8];    % 深蓝（原点）
color_1 = [0.8, 0.4, 0.2];    % 暖红
color_2 = [0.2, 0.7, 0.5];    % 浅绿
color_3 = [0.6, 0.4, 0.8];    % 淡紫
color_A = [0.3, 0.3, 0.3];    % 深灰

% 节点样式（填充+边缘+立体）
plot(S(1), S(2), 'o', ...
    'MarkerSize', 12, 'MarkerFaceColor', color_S, 'MarkerEdgeColor', 'white', ...
    'LineWidth', 1.5, 'DisplayName', 'S（原点）');

plot(node1(1), node1(2), 'o', ...
    'MarkerSize', 12, 'MarkerFaceColor', color_1, 'MarkerEdgeColor', 'white', ...
    'LineWidth', 1.5, 'DisplayName', '节点1');

plot(node2(1), node2(2), 'o', ...
    'MarkerSize', 12, 'MarkerFaceColor', color_2, 'MarkerEdgeColor', 'white', ...
    'LineWidth', 1.5, 'DisplayName', '节点2');

plot(node3(1), node3(2), 'o', ...
    'MarkerSize', 12, 'MarkerFaceColor', color_3, 'MarkerEdgeColor', 'white', ...
    'LineWidth', 1.5, 'DisplayName', '节点3');

plot(A(1), A(2), 'o', ...
    'MarkerSize', 12, 'MarkerFaceColor', color_A, 'MarkerEdgeColor', 'white', ...
    'LineWidth', 1.5, 'DisplayName', '节点A');

% ===================== 4. 绘制连线（分层+粗细区分） =====================
% 主连线（S与各节点，加粗）
plot([S(1), node1(1)], [S(2), node1(2)], '-', 'Color', color_S, 'LineWidth', 2);
plot([S(1), node2(1)], [S(2), node2(2)], '-', 'Color', color_S, 'LineWidth', 2);
plot([S(1), node3(1)], [S(2), node3(2)], '-', 'Color', color_S, 'LineWidth', 2);

% 次连线（其他节点间，稍细）
plot([node1(1), node2(1)], [node1(2), node2(2)], '-', 'Color', color_1, 'LineWidth', 1.2);
plot([node1(1), A(1)], [node1(2), A(2)], '-', 'Color', color_1, 'LineWidth', 1.2);
plot([node2(1), A(1)], [node2(2), A(2)], '-', 'Color', color_2, 'LineWidth', 1.2);
plot([node3(1), node2(1)], [node3(2), node2(2)], '-', 'Color', color_3, 'LineWidth', 1.2);
plot([node3(1), A(1)], [node3(2), A(2)], '-', 'Color', color_3, 'LineWidth', 1.2);

% ===================== 5. 标注距离（移除Padding，低版本兼容） =====================
% 标注样式统一配置（删除Padding，保留其他美化属性）
textOpt = struct(...
    'FontName', '微软雅黑', ...
    'FontSize', 10, ...
    'FontWeight', 'bold', ...
    'BackgroundColor', [1,1,1], ... % 白色背景框
    'EdgeColor', [0.8,0.8,0.8]);  % 浅灰边框

% S相关标注
text(mean([S(1), node1(1)])+5, mean([S(2), node1(2)])+5, '160', textOpt);
text(mean([S(1), node2(1)]), mean([S(2), node2(2)])-8, '140', textOpt);
text(mean([S(1), node3(1)])+5, mean([S(2), node3(2)])-10, '160', textOpt);

% 其他标注
text(mean([node1(1), node2(1)])+5, mean([node1(2), node2(2)]), '130', textOpt);
text(mean([node1(1), A(1)])-8, mean([node1(2), A(2)])+5, '160', textOpt);
text(mean([node2(1), A(1)])+8, mean([node2(2), A(2)]), '140', textOpt);
text(mean([node3(1), node2(1)])+5, mean([node3(2), node2(2)]), '130', textOpt);
text(mean([node3(1), A(1)])-8, mean([node3(2), A(2)])-5, '160', textOpt);

% ===================== 6. 图例美化 =====================
lgd = legend('Location', 'bestoutside');
lgd.FontName = '微软雅黑';
lgd.FontSize = 10;
lgd.EdgeColor = 'none'; % 去掉图例边框
% 低版本Matlab无CellSpacing，注释掉该属性
% lgd.CellSpacing = 5;    

% ===================== 7. 节点名称标注（精准+美观） =====================
nameOpt = struct(...
    'FontName', '微软雅黑', ...
    'FontSize', 11, ...
    'FontWeight', 'bold', ...
    'Color', [0.1,0.1,0.1]);

text(S(1)-8, S(2)+8, 'S', nameOpt);
text(node1(1)+8, node1(2)+8, '1', nameOpt);
text(node2(1)+8, node2(2)+8, '2', nameOpt);
text(node3(1)+8, node3(2)-8, '3', nameOpt);
text(A(1)+8, A(2)+8, 'A', nameOpt);

hold off;

% 可选：导出高清图（300dpi）
% print(gcf, '地理位置关系图', '-dpng', '-r300');