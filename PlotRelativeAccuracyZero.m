
HA=[[HistoriedRelativeAccuracy.allcat]; [HistoriedRelativeAccuracy.cat]; ...
    [HistoriedRelativeAccuracy.dom]; [HistoriedRelativeAccuracy.etd]; ...
    [HistoriedRelativeAccuracy.ed]; [HistoriedRelativeAccuracy.morph]; ...
    [HistoriedRelativeAccuracy.num]; [HistoriedRelativeAccuracy.tec]; ...
    [HistoriedRelativeAccuracy.last]; ...
    [HistoriedRelativeAccuracy.median]; ...
    [HistoriedRelativeAccuracy.mode]; ...
    [HistoriedRelativeAccuracy.min]; ...
    [HistoriedRelativeAccuracy.max]];


h=bar(HA);
%l=legend([Result(:,1,:,1).VolcanoDatasetName]);
labels=[{'VEI $\ge$ 0'},{'VEI $\ge$ 1'},{'VEI $\ge$ 2'},{'VEI $\ge$ 3'}, ...
    {'VEI $\ne$ 2'}];
l=legend(labels,'interpreter','latex');
%x=[{'All Categories'},{'Tectonic Setting'},{'Dominant Rock Type'},{'Morphology'}];
Attributes = {'All Attributes','Categorical Data','Dominant Rock Type', ...
    'Repose Time','Eruption Duration','Morphology','All Numerical Data','Tectonic Setting' ...
    'Last VEI','Median VEI','Mode VEI','Min VEI','Max VEI'};
set(gca,'XTickLabel',Attributes)

% errorbars
hold on
for ih=1:length(h),
        errorbar(h(ih).XEndPoints,h(ih).YData, ...
            HistoriedErr(ih)+0*h(ih).YData,'Color',[0.7 0.7 0.7], ...
            'LineWidth',2,'CapSize',4,'LineStyle','none')
end
hold off

% population numbers
for ih=1:length(h),
        I=(h(ih).YData)>0;
      % positive valuaes aligned to top of bars
        text(h(ih).XEndPoints(I),h(ih).YData(I),num2cell(HN(I,ih)), ...
        'HorizontalAlignment','right','VerticalAlignment','middle', ...
        'Color', 'w', 'Rotation', 90,'Fontsize',12)
% negative values aligned to bottom of bars
        text(h(ih).XEndPoints(~I),h(ih).YData(~I),num2cell(HN(~I,ih)), ...
        'HorizontalAlignment','left','VerticalAlignment','middle', ...
        'Color', 'w', 'Rotation', 90,'Fontsize',12)
 
end
% separate simple predictors
line(([1 1]*h(5).XEndPoints(8)+h(1).XEndPoints(9))/2,[0 1], ...
    'Color','k','LineStyle','--')
ylim(RelAccYlims)
legend(l.String{1:5})
ylabel('Accuracy Gain')