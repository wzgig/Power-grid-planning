%% =========================
%  方案一（最终美观版）
%  固定顺序：发电厂A -> 变电所1 -> 变电所2 -> 变电所3
% =========================

%% 站点名称（按你要求固定顺序）
stations = {'发电厂A','变电所1','变电所2','变电所3'};

%% 数据（题目给定 I/II；III = 100 - I - II）
I  = [35 10 30  0];
II = [25 35 30 25];
III = 100 - I - II;

X = [I(:), II(:), III(:)];
typeNames = {'I类','II类','III类'};

%% ——可调参数——
labelThreshold = 7;            % 段内标注阈值：小于该值不在段内写，避免拥挤
showRightNoteForSmall = true;  % 小段标签是否统一放到右侧说明

%% 配色（稳重、报告友好）
C = [  0   114  178;   % 蓝
      213   94    0;   % 橙
       86  180  120] / 255; % 绿

%% 画图
figure('Color','w','Position',[80 80 1100 540]);  % 宽一点，右侧要放小段说明

b = barh(X,'stacked','BarWidth',0.62,'LineWidth',0.8);
for k = 1:numel(b)
    b(k).FaceColor = C(k,:);
    b(k).EdgeColor = [1 1 1];  % 白色分隔线（关键：层次清晰）
    b(k).LineWidth = 1.1;
end

ax = gca;
ax.FontName = 'Microsoft YaHei';
ax.FontSize = 12;
ax.LineWidth = 1.0;
ax.Box = 'off';
ax.TickDir = 'out';
ax.XLim = [0 100];
ax.XTick = 0:10:100;

% y轴：按你要求从上到下显示 A,1,2,3
ax.YTick = 1:numel(stations);
ax.YTickLabel = stations;
set(ax,'YDir','reverse');      % 让第1个出现在最上面（barh默认第1个在下）

ax.XGrid = 'on';
ax.GridAlpha = 0.12;
ax.YGrid = 'off';

% 参考线（可选，增强读图）
xline(50,'--','LineWidth',0.8,'Alpha',0.25);

xlabel('占比（%）','FontWeight','bold');
title('各站点负荷类型占比（100%堆叠条形图）','FontWeight','bold');

% 图例：置顶横排
lgd = legend(typeNames,'Location','northoutside','Orientation','horizontal');
lgd.Box = 'off';

%% 左侧小注释（加码）
noteLeft = sprintf('注：III类占比 = 100%% − I类 − II类（由题目给定I、II类计算）');
text(-2, 0.2, noteLeft, ...
    'Units','data', 'HorizontalAlignment','left','VerticalAlignment','top', ...
    'FontName','Microsoft YaHei','FontSize',10,'Color',[0.35 0.35 0.35]);

%% 段内标注：自动黑/白字（根据背景亮度）
luma = 0.2126*C(:,1) + 0.7152*C(:,2) + 0.0722*C(:,3); % 亮度估计

hold on;
cumX = cumsum(X,2);

for i = 1:size(X,1)
    for j = 1:size(X,2)
        val = X(i,j);
        if val >= labelThreshold
            xCenter = cumX(i,j) - val/2;

            % 亮背景用黑字，暗背景用白字
            txtColor = [1 1 1];
            if luma(j) > 0.62, txtColor = [0 0 0]; end

            text(xCenter, i, sprintf('%d%%', round(val)), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontName','Microsoft YaHei', ...
                'FontSize',11, ...
                'FontWeight','bold', ...
                'Color',txtColor);
        end
    end
end

%% 小段标签：右侧统一说明（加码，更整洁）
if showRightNoteForSmall
    xRight = 102.0;   % 右侧说明起始位置
    for i = 1:size(X,1)
        smallParts = find(X(i,:) < labelThreshold & X(i,:) > 0);
        if ~isempty(smallParts)
            pieces = strings(1,numel(smallParts));
            for t = 1:numel(smallParts)
                jj = smallParts(t);
                pieces(t) = sprintf('%s %d%%', typeNames{jj}, round(X(i,jj)));
            end
            text(xRight, i, strjoin(pieces,'  |  '), ...
                'HorizontalAlignment','left', ...
                'VerticalAlignment','middle', ...
                'FontName','Microsoft YaHei', ...
                'FontSize',10, ...
                'Color',[0.25 0.25 0.25]);
        end
    end
    ax.XLim = [0 120];     % 给右侧说明留白
    ax.XTick = 0:10:100;   % 刻度仍显示到 100
end

hold off;

%% 右下角数据来源脚注（加码）
% 用 annotation 更稳，不随坐标缩放
foot = '数据来源：课程设计题目表1-1（I类、II类）；III类=100%-I类-II类（计算）。';
annotation('textbox',[0.02 0.01 0.96 0.05], ...
    'String',foot, 'EdgeColor','none', ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
    'FontName','Microsoft YaHei','FontSize',10,'Color',[0.35 0.35 0.35]);

%% 导出（可选）
% exportgraphics(gcf,'负荷类型占比_方案一_固定顺序.png','Resolution',300);
