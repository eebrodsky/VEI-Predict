function MakePieChart(Result,Model, attribute)
% creates pie charts comparing the distribution of attribute in 
% populations correctly predicted by modeltype versus full dataset
% Result structure
FS= 18;
R=Result;
s=['R.AttributeTable.Tbl_VEI_allcat.' attribute];
eval(['fulldataset=' s ';']);
tbl=R.AttributeTable.Tbl_VEI_allcat;
trainpredict=Model.predictFcn(tbl);
Icorrect=(trainpredict==tbl.currentVEI);
correct=fulldataset(Icorrect);

[full_counts,full_labels]=hist(categorical(fulldataset));
[correct_counts]=hist(categorical(correct),full_labels);

subplot(1,2,1)
pie(full_counts);
legend(full_labels,'location','best','FontSize',FS);
title('Full Dataset')
subplot(1,2,2)
pie(correct_counts);
(title('Correctly Predicted'))
%legend(correct_labels,'location','best','FontSize',FS);