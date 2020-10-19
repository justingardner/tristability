function [] = Figure1_ExperimentalDesign()
%This function makes a diagram plotting how each trial type proceeds
%through orientation disparities over time.
mSize = 7;
LTM = .2;
LT = .2;
fontSize = 10;
titleFontSize = 12;

h = figure('rend','painter','pos',[10 10 900 1650]);
SlowOutColor = [0.5686    0.2353    0.3333];%Dark Red
SlowInColor = [0.2157    0.3373    0.4549];%Dark Blue
FastOutColor = [0.8235    0.4392    0.5529];%Light Red
FastInColor = [0.5412    0.6392    0.7647];%Light Blue

plot([10 40],[160 140],'-','color',FastOutColor,'LineWidth', LT,'MarkerEdgeColor','w');
hold on;
plot([40 10],[160 140],'-','color',FastInColor,'LineWidth', LT,'MarkerEdgeColor','w');
plot([10 40],[160 0],'-','color',SlowOutColor,'LineWidth', LT,'MarkerEdgeColor','w');
plot([40 10],[160 0],'-','color',SlowInColor,'LineWidth', LT,'MarkerEdgeColor','w');


legendText = {' Fast Out',' Fast In',' Slow Out',' Slow In'};
legendColors = {FastOutColor,FastInColor,SlowOutColor,SlowInColor};
legendX = 30;
legendY = 10;
legendSpacing = -5;
legendGap = 2;
legendLineSegSize = 10;
for i=1:4
    plot([legendX-legendLineSegSize-legendGap legendX-legendGap],[legendY-(legendSpacing*i) legendY-(legendSpacing*i)],'lineWidth',LT,'color',legendColors{i});
    text(legendX,legendY-(legendSpacing*i),legendText{i},'VerticalAlignment','middle','HorizontalAlignment','left','FontAngle','oblique','FontSize',fontSize);
end

line([10 40],[160 140],'Color',FastOutColor,'MarkerEdgeColor',[.999 .999 .999],'MarkerSize', mSize,'LineWidth', LTM);
line([40 10],[160 140],'Color',FastInColor,'MarkerEdgeColor',[.999 .999 .999],'MarkerSize', mSize,'LineWidth', LTM);
line([10 40],[160 0],'Color',SlowOutColor,'MarkerEdgeColor',[.999 .999 .999],'MarkerSize', mSize,'LineWidth', LTM);
line([40 10],[160 0],'Color',SlowInColor,'MarkerEdgeColor',[.999 .999 .999],'MarkerSize', mSize,'LineWidth', LTM);



%xaxis(20,40);
xlabel('Angle between Gabors (degrees)');
ylabel('Trial time (seconds)');
drawPublishAxis('xTickOffset',-.005,'yTickOffset',.5,'whichAxis','both','lineWidth',1,'labelFontSize',fontSize,'titleStr','b','xAxisMajorTickLen',-1/64,'xAxisMinorTickLen',-1/128,'xAxisOffset',1/10,'xAxisMin',0,'xAxisMax',50,'xTick',0:10:50,'xTickLabel',0:10:50,'xTickLabelSigfigs',0,'yTickLabelSigfigs',0,'yAxisMin',0,'yAxisMax',160,'yAxisMajorTickLen',-1/64,'yTick',0:20:160,'yTickLabel',160:-20:0)
xlim([0 50]);

set(h,'PaperUnits','inches');
x_width = 3.5;
y_width = 10;
set(h,'PaperPosition',[.5 .5 x_width y_width]);
saveas(h,'Figure1_ExperimentalDesign_raw.pdf')


end

