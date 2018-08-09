mymap = ones(1001,3);
mymap(1:501,1) = [0:0.002:1]';
mymap(1:501,2) = mymap(1:501,1);
mymap(501:1001,2) = flip(mymap(1:501,1));
mymap(501:1001,3) = mymap(501:1001,2);
% mymap(491:511,:) = 1;
mymap = mymap.^0.5;
F = {'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','O1','O2','O3','O4','O5','O6','O7','O8','O9','O10','O11'};

fig1 = openfig('G by AFA_ALDM.fig');
ax = gca;
title('');
colormap(mymap)
set(fig1, 'Position', [0 0 1920 1080]);
colorbar('westoutside')
caxis([-2.5 2.5]);
im.AlphaData = 0.5;
% xlabel('Feature');
% xticks([1:1:G_n]);
xlabel('');
xticklabels(F);
set(gca, 'XTickLabel', F);
ylabel('Index');
% ax.FontSize = 20;

fig2 = openfig('G by AFA_ALDM_30.fig');
% left_color = [.5 .5 0];
% right_color = [0 .5 .5];
% set(fig,'ColorOrder',[left_color; right_color]);
title('');
set(fig2, 'Position', [0 0 1920 1080]);
colormap(mymap)
colorbar('westoutside')
caxis([-2.5 2.5])
im.AlphaData = .5;
% xlabel('Feature');
xticks([1:1:G_n]);
xticklabels(F);

yyaxis left
ax = gca;
ax.YDir = 'reverse';
ylabel('Index','Color','k');
yticks([1:1:30]);
yticklabels([1:1:30]);
% yyaxis(ax,'right','color','k');
yyaxis right
ax.YAxis(2,1).Color = 'black';
ax.YDir = 'reverse';
% ylabel('Chromosome and Location','Color','k');
ylabel('P Value','Color','k');
% ax.YAxis(2,1).yticks([1:1:30]);
ax.YAxis(2,1).TickLabels = G_sto(1:30,25);
ax.YAxis(2,1).Limits = ax.YAxis(1,1).Limits;
ax.YAxis(2,1).TickValues = ax.YAxis(1,1).TickValues;
ax.FontSize = 20;
% ax.YAxis(3,1) = ax.YAxis(2,1);