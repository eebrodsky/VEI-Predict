h=bar([[UnHistoriedRelativeAccuracy.cat]; [UnHistoriedRelativeAccuracy.tec]; [UnHistoriedRelativeAccuracy.dom]; [UnHistoriedRelativeAccuracy.morph]]);
x=[{'All Categories'},{'Tectonic Setting'},{'Dominant Rock Type'},{'Morphology'}];
set(gca,'XTickLabel',x)
ylabel('Accuracy Gain')

% errorbars
if (errorbarflag)
hold on
for ih=1:length(h),
        errorbar(h(ih).XEndPoints,h(ih).YData, ...
            HistoriedErr(ih)+0*h(ih).YData,'Color',[0.7 0.7 0.7], ...
            'LineWidth',2,'CapSize',4,'LineStyle','none')
end
hold off
end

for ih=1:length(h),
    text(h(ih).XEndPoints,h(ih).YData,num2cell(HN(:,ih)), ...
        'HorizontalAlignment','right','VerticalAlignment','middle', ...
        'Color', 'w', 'Rotation', 90,'Fontsize',11)
end
ylim(RelAccYlimsPos)
legend('off')
text(0.01,0.9,'\bf b','Units','normalized','FontSize',24)

