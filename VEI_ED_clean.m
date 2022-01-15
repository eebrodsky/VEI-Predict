%% Cleaned code for Predicting VEI from Global Database 
clear


%%%%%% Download Data %%%%%%%%

%Download the eruption data from the data file
eruption_data =readtable('GVP_Eruption_Results_Clean_3.csv');

%Download the volcano data catelog from the data file
volcano_data =readtable('New_GVP_Volcano_List_Holocene.csv');

%% Creating a Baseline for Accuracy and RMSE 
%%%%%%% Initiation Section %%%%%%%%
Nbstrap=500; % Number of bootstraps for baseline calculations
% naming conventions will be: VolcanoDatasetName=[StartTime{itime} '_' VEIType{i_bol} '_VEI' iVEI ];
StartTime=[{'1500'},{'all'}];
VEIType=[{'threshold','remove'}];
NResult=4*2*2*2; % Maximum number of Result records

%Varying VEI thresholds  
for iVEI=0:3
    %Varying Time thresholds
    for itime = 1:2
        
        if itime == 1 
            eruption_data_time = 1500;
        else
            eruption_data_time = min(eruption_data.StartYear);
        end 
        %Changing between VEI thresholds and Removing specific VEI (only
        %removing VEI =2
        if (iVEI==2) 
            Nbol=2;
        else
            Nbol=1;
        end
            for i_bol = 1:Nbol
            
            if i_bol == 1 
                eruption_data_complete = eruption_data(eruption_data.StartYear >= eruption_data_time,:);
                eruption_data_clean = eruption_data_complete(eruption_data_complete.VEI >= iVEI,:);      
            else 
                eruption_data_complete = eruption_data(eruption_data.StartYear >= eruption_data_time,:);
                eruption_data_clean = eruption_data_complete(eruption_data_complete.VEI ~= iVEI,:);
            end
            % set up the name of this dataset
            VolcanoDatasetName=cellstr([StartTime{itime} ' ' VEIType{i_bol} ' VEI' num2str(iVEI) ]);

            %Remove NAN values from Eruption_data
            eruption_data_clean_ind = find(isnan(eruption_data_clean.VEI) == 0);
            eruption_data_clean_nan = eruption_data_clean(eruption_data_clean_ind,:);

            %find unique Volcano names and numbers  
            [VolcanoList,~,IVolcanoList]=unique(eruption_data_clean_nan.VolcanoName);
            [VolcanoList_cat,~,IVolcanoList_cat]=unique(volcano_data.VolcanoNumber);

            %total number of unique volcano names.
            N=length(VolcanoList);
            Volcano=[];
            %Looping to File each eruption with a unique Volcano
            for i=1:N

                %finds all of the indices of each volcano from the eruption_data set 
                J=find(IVolcanoList==i);
                
                %Creates a structure of called Volcano and runs through all of the
                %unique names of volcanos and creates a row in a column called name in 
                %the volcano structure and it draws these name from the object
                %volcanoList that was created above from the erutpion _data frame
                Volcano(i).name=VolcanoList(i);
                
                         
                
                %Loop through all the eruptions corresponding to a given
                %volcano
                for j=1:length(J),

                    if j==1
                        
                        %This creates a new column in the structure volcano called
                        %number and draws that number directly from the eruption data
                        %frame if the j index is 1 and then assigns it to the unqiue
                        %volcano name.
                        Volcano(i).number=eruption_data_clean_nan(J(j),:).VolcanoNumber;

                    end
                    
                    %Row goes into the eruption data frame and takes all of the
                    %eruption events index by J for each loop of j which is the length
                    %of the number of indexes that each volcano has for the number of
                    %eruptions that have occured in the holocene.
                    Row=eruption_data_clean_nan(J(j),:);
                    
                    %Creates a column that in the structure volcano called eruption and
                    %column is made up of structures that contain all the eruption 
                    %events for the unique volcano names. This uses the function
                    %created below so that the structure collects the following
                    %attributes from the eruption data set.
                    Volcano(i).eruption(j)=FileEruption(Row);
                end
            end
% resort  eruptions in descending year order with
% the most recent (called "current") eruption first
            for k=1:length(Volcano),
                s=[Volcano(k).eruption.startyear];
                [~,I]=sort(s,'descend');
                resorted=Volcano(k).eruption(I);
                Volcano(k).eruption=resorted;
            end
            
            %CheckVolcanoYearOrder;
  % pull out attributes from volcano_data
             for k=1:length(Volcano),
                 L=find(Volcano(k).number==volcano_data.VolcanoNumber);
                 if (length(L)>1)
                     ['Warning: non-unique volcano identification']
                 end
                 Volcano(k).rocktype=volcano_data(L,:).DominantRockType;
                 Volcano(k).tectonicsetting=volcano_data(L,:).TectonicSetting;
                 Volcano(k).primaryvolcanotype=volcano_data(L,:).PrimaryVolcanoType;
                 Volcano(k).latitude=volcano_data(L,:).DominantRockType;
                 Volcano(k).longitude=volcano_data(L,:).TectonicSetting;
             end
                    
            
            
            %%%%%%%%%%Creating a Baseline for Historied Volcanoes%%%%%%%%%% 
            
            for k=1:length(Volcano)
                
                % Find Historied Volcanoes
                if (length(Volcano(k).eruption)>1)

                    Volcano(k).currentVEI=Volcano(k).eruption(1).vei;
                    
                      %Calculate previous (for persistence), median, mode, min, max 
                    Volcano(k).SimplePredictors.last = Volcano(k).eruption(2).vei;
                    Volcano(k).SimplePredictors.median = round(median([Volcano(k).eruption(2:end).vei]));
                    Volcano(k).SimplePredictors.mode = (mode([Volcano(k).eruption(2:end).vei]));
                    Volcano(k).SimplePredictors.min = (min([Volcano(k).eruption(2:end).vei]));
                    Volcano(k).SimplePredictors.max = (max([Volcano(k).eruption(2:end).vei]));
                    
                else
                    %This makes all of the volcanoes with 1 eruption labeled NaN because
                    %we would not be able to have a test group and then a value to predict
                    Volcano(k).currentVEI=NaN;
                    
                end
            end
            
            currentVEI=[Volcano.currentVEI];
            Ihistory=isfinite(currentVEI);
            Iunhistory=isnan(currentVEI);
            VolcanoHistory=Volcano(Ihistory);
            VolcanoUnhistory=Volcano(Iunhistory);
            
% process and calculate baseline accuracy for Historied Volcanos 
            ihist=1;
            Result(iVEI+1,itime,i_bol,ihist).Volcano=VolcanoHistory;
            Result(iVEI+1,itime,i_bol,ihist).VolcanoDatasetName=VolcanoDatasetName;
            Result(iVEI+1,itime,i_bol,ihist).HistoryFlag=1;
            Result(iVEI+1,itime,i_bol,ihist).eruption_data=eruption_data_clean_nan;

            currentVEI=[VolcanoHistory.currentVEI];
            Nvolc=length(currentVEI);
            for j_i=1:Nbstrap
                I=randi(Nvolc,1,Nvolc);
                Synth=currentVEI(I);
                Accuracy_boot(j_i)=(sum(Synth==currentVEI))/Nvolc;
                RMSE_boot(j_i) = (sum(((Synth-currentVEI).^2)/Nvolc)).^(1/2);

            end

            % 0-rule accuracy using minimum VEI in dataset
            ZeroRule=mode(currentVEI);
            Result(iVEI+1,itime,i_bol,ihist).ZeroRule=mean(currentVEI==ZeroRule);

            %Accuracy Baseline based on the distribution 
            Result(iVEI+1,itime,i_bol,ihist).BaselineAccuracy= mean(Accuracy_boot); 

            %Accuracy Baseline Standard Deviation 
            Result(iVEI+1,itime,i_bol,ihist).BaselineStd = std(Accuracy_boot);
            
            %RMSE Baseline based on the distribution
            Result(iVEI+1,itime,i_bol,ihist).BaselineRMSE = mean(RMSE_boot);

            %%%%%%%%%Process and create a Baseline for Unhistoried Volcanoes%%%%%%%%%
            ihist=2;
            Result(iVEI+1,itime,i_bol,ihist).Volcano=VolcanoUnhistory;
            Result(iVEI+1,itime,i_bol,ihist).VolcanoDatasetName=VolcanoDatasetName;
            Result(iVEI+1,itime,i_bol,ihist).HistoryFlag=0;

            VEIlist=[];
            for i=1:length(VolcanoUnhistory),
                VEIlist(i)=VolcanoUnhistory(i).eruption.vei;
            end
            Nvolc=length(VEIlist);
            Result(iVEI+1,itime,i_bol,ihist).currentVEI=VEIlist;

            
            %Calculating the Accuracy baseline
            for j_i=1:Nbstrap
                I=randi(Nvolc,1,Nvolc);
                Synth=VEIlist(I);
                RMSE_boot_unhist(j_i) = (sum(((Synth-VEIlist).^2)/Nvolc)).^(1/2);

                Accuracy_boot_unhist(j_i)=(sum(Synth==VEIlist))/Nvolc;
            end

            %Accuracy Baseline based on the distribution 
            Result(iVEI+1,itime,i_bol,ihist).BaselineAccuracy= mean(Accuracy_boot_unhist); 

            % 0-rule accuracy using minimum VEI in dataset
            ZeroRule=mode(VEIlist);
            Result(iVEI+1,itime,i_bol,ihist).ZeroRule=mean(VEIlist==ZeroRule);

            %Accuracy Baseline Standard Deviation 
            Result(iVEI+1,itime,i_bol,ihist).BaselineStd = std(Accuracy_boot_unhist);
            
            %RMSE Baseline based on the distribution
            Result(iVEI+1,itime,i_bol,ihist).BaselineRMSE = mean(RMSE_boot_unhist);

            end
        
    end 
end
% save intermediate results
filename='Result';
save(filename,'Result');

%%
% Do simple predictors on historied datasets (HistoryFlag=1)
% look through each of the simple predictors and store the accuracy in the
% Result structure
for i=1:NResult % looping linearly through the 
    %Result array to avoid excessive for loop indices
    if(Result(i).HistoryFlag) % only do this for historied volcanoes
        PredictorNames=fieldnames(Result(i).Volcano(1).SimplePredictors);
        Predictors=[Result(i).Volcano.SimplePredictors];
        for j=1:length(PredictorNames),
            P=strcat('[Predictors.', PredictorNames(j), ']');
            VEIPredict=eval(P{1});
            CurrentVEI=[Result(i).Volcano.currentVEI];
            Accuracy=mean(CurrentVEI==VEIPredict);
            RMSE=(sum(((VEIPredict-CurrentVEI).^2)/length(CurrentVEI))).^(1/2);
            StoreStr = strcat('Result(i).Accuracy.', ...
                PredictorNames(j), '= Accuracy;');
            eval(StoreStr{1})
            StoreStr = strcat('Result(i).RelativeAccuracy.', ...
                PredictorNames(j), '= Accuracy-Result(i).BaselineAccuracy;');
            eval(StoreStr{1})
            StoreStr = strcat('Result(i).RelativeAccuracyZero.', ...
                PredictorNames(j), '= Accuracy-Result(i).ZeroRule;');
            eval(StoreStr{1})
            StoreStr = strcat('Result(i).RMSE.', ...
                PredictorNames(j), '= RMSE;');
            eval(StoreStr{1})
            Result(i).currentVEI=CurrentVEI;
            Result(i).NumberVolcano.Simple=length(CurrentVEI);
        end
    end
end
% save intermediate results
filename='Result';
save(filename,'Result');

%% Creating Tables of Attributes 

% VEI Thresholds 
for iResult=1:NResult,
    if (~isempty(Result(iResult).Volcano))
        [Tbl_VEI_allcat] = Create_table_attributes_ED_clean(Result(iResult).Volcano, Result(iResult).HistoryFlag);
        Result(iResult).AttributeTable.Tbl_VEI_allcat=Tbl_VEI_allcat;
    end
end
%%
retrainflag = true;

if retrainflag == 1
    tic
  
    for iResult = 1:NResult,
        if (Result(iResult).HistoryFlag)&(~isempty(Result(iResult).Volcano))
            [iResult toc]
         table_names_mach_learn = ["allcat","cat","dom","etd","ed","morph","num","tec"];
        Tbl_VEI_allcat=Result(iResult).AttributeTable.Tbl_VEI_allcat;
        %Historied Volcano Model and Table of Attributes 
        for k = 1:length(table_names_mach_learn)
            if k == 1
                % Everything called: allcat
                    tbl=Tbl_VEI_allcat;
                    Imissing=cellfun(@(x) isequal(x,'999'),tbl.RockType);
                    tbl=tbl(~Imissing,:);
                    Imissing=isnan(tbl.eruptiontimediff);
                    tbl=tbl(~Imissing,:);
                    Imissing=isnan(tbl.eruptionduration);
                    tbl=tbl(~Imissing,:);

                    p=tbl.Properties.VariableNames;
                  
                    I=cellfun(@(x) isequal(x,'currentVEI')|isequal(x,'name'),p);
                    PredictVect=p(~I);
                    [model,accuracy] = trainClassifier_Generic(tbl,PredictVect);
                    %f(k)=parfeval(@trainClassifier_Generic,2,tbl,PredictVect);
                    %[model,accuracy]=fetchOutputs(f(k));

                    trainpredict=model.predictFcn(tbl);
                    trainaccuracy=mean(trainpredict==tbl.currentVEI);
                   
                    Result(iResult).NumberVolcano.allcat=height(tbl);
                    Result(iResult).Model.allcat=model;
                    Result(iResult).Accuracy.allcat=accuracy;
                    Result(iResult).RelativeAccuracy.allcat=accuracy-Result(iResult).BaselineAccuracy;
                    Result(iResult).RelativeAccuracyZero.allcat=accuracy-Result(iResult).ZeroRule;
                    Result(iResult).TrainAccuracy.allcat=trainaccuracy;
                    Result(iResult).TrainRelativeAccuracy.allcat=trainaccuracy-Result(iResult).BaselineAccuracy;

            elseif k == 2
                % all categorical called: cat
                    tbl=Tbl_VEI_allcat;
                    Imissing=cellfun(@(x) isequal(x,'999'),tbl.RockType);
                    tbl=tbl(~Imissing,:);

                    p=tbl.Properties.VariableNames;
                    I=cellfun(@(x) isequal(x,'RockType')|isequal(x,'TectonicSetting')|isequal(x,'morphology'),p);
                    PredictVect=p(I);
                    [model,accuracy] = trainClassifier_Generic(tbl, PredictVect);
             %       f(k)=parfeval(@trainClassifier_Generic,2,tbl,PredictVect);
              %      [model,accuracy]=fetchOutputs(f(k));

                    trainpredict=model.predictFcn(tbl);
                    trainaccuracy=mean(trainpredict==tbl.currentVEI);

                 
                Result(iResult).NumberVolcano.cat=height(tbl);
                Result(iResult).Model.cat=model;
                Result(iResult).Accuracy.cat=accuracy;
                Result(iResult).RelativeAccuracy.cat=accuracy-Result(iResult).BaselineAccuracy;
                Result(iResult).RelativeAccuracyZero.cat=accuracy-Result(iResult).ZeroRule;
                Result(iResult).TrainAccuracy.cat=trainaccuracy;
                Result(iResult).TrainRelativeAccuracy.cat=trainaccuracy-Result(iResult).BaselineAccuracy;


            elseif k == 3
                % Rocktype called: dom
                   tbl=Tbl_VEI_allcat;
                    Imissing=cellfun(@(x) isequal(x,'999'),tbl.RockType);
                    tbl=tbl(~Imissing,:);

                    p=tbl.Properties.VariableNames;
                    I=cellfun(@(x) isequal(x,'RockType'),p);
                    
                    PredictVect=p(I);
                    [model,accuracy] = trainClassifier_Generic(tbl, PredictVect);
                %    f(k)=parfeval(@trainClassifier_Generic,2,tbl,PredictVect);
                %    [model,accuracy]=fetchOutputs(f(k));
                    trainpredict=model.predictFcn(tbl);
                    trainaccuracy=mean(trainpredict==tbl.currentVEI);
                 
                Result(iResult).NumberVolcano.dom=height(tbl);
                Result(iResult).Model.dom=model;
                Result(iResult).Accuracy.dom=accuracy;
                Result(iResult).RelativeAccuracy.dom=accuracy-Result(iResult).BaselineAccuracy;
                Result(iResult).RelativeAccuracyZero.dom=accuracy-Result(iResult).ZeroRule;
                Result(iResult).TrainAccuracy.dom=trainaccuracy;
                Result(iResult).TrainRelativeAccuracy.dom=trainaccuracy-Result(iResult).BaselineAccuracy;


            elseif k == 4
                % Eruption Time difference   called: etd
                    tbl=Tbl_VEI_allcat;
                      Imissing=isnan(tbl.eruptiontimediff);
                    tbl=tbl(~Imissing,:);
                    p=tbl.Properties.VariableNames;
                    I=cellfun(@(x) isequal(x,'eruptiontimediff'),p);
                    PredictVect=p(I);
                    [model,accuracy] = trainClassifier_Generic(tbl, PredictVect);
                   %  f(k)=parfeval(@trainClassifier_Generic,2,tbl,PredictVect);
                   % [model,accuracy]=fetchOutputs(f(k));
                    trainpredict=model.predictFcn(tbl);
                    trainaccuracy=mean(trainpredict==tbl.currentVEI);

                    Result(iResult).NumberVolcano.etd=height(tbl);
                    Result(iResult).Model.etd=model;
                    Result(iResult).Accuracy.etd=accuracy;
                    Result(iResult).RelativeAccuracy.etd=accuracy-Result(iResult).BaselineAccuracy;
                    Result(iResult).RelativeAccuracyZero.etd=accuracy-Result(iResult).ZeroRule;
                    Result(iResult).TrainAccuracy.etd=trainaccuracy;
                    Result(iResult).TrainRelativeAccuracy.etd=trainaccuracy-Result(iResult).BaselineAccuracy;

            elseif k == 5
                % Duration called: ed
                    tbl=Tbl_VEI_allcat;
                       Imissing=isnan(tbl.eruptionduration);
                    tbl=tbl(~Imissing,:);
                    p=tbl.Properties.VariableNames;
                    I=cellfun(@(x) isequal(x,'eruptionduration'),p);
                    PredictVect=p(I);
                   
                    [model,accuracy] = trainClassifier_Generic(tbl, PredictVect);
                   % f(k)=parfeval(@trainClassifier_Generic,2,tbl,PredictVect);
                   % [model,accuracy]=fetchOutputs(f(k));
                    trainpredict=model.predictFcn(tbl);
                    trainaccuracy=mean(trainpredict==tbl.currentVEI);


                    Result(iResult).NumberVolcano.ed=height(tbl);
                    Result(iResult).Model.ed=model;
                    Result(iResult).Accuracy.ed=accuracy;
                    Result(iResult).RelativeAccuracy.ed=accuracy-Result(iResult).BaselineAccuracy;
                    Result(iResult).RelativeAccuracyZero.ed=accuracy-Result(iResult).ZeroRule;
                    Result(iResult).TrainAccuracy.ed=trainaccuracy;
                    Result(iResult).TrainRelativeAccuracy.ed=trainaccuracy-Result(iResult).BaselineAccuracy;


            elseif k == 6
                % morphology called: morph
                    tbl=Tbl_VEI_allcat;
                    p=tbl.Properties.VariableNames;
                    I=cellfun(@(x) isequal(x,'morphology'),p);
                    PredictVect=p(I);
                    [model,accuracy] = trainClassifier_Generic(tbl, PredictVect);
                   % f(k)=parfeval(@trainClassifier_Generic,2,tbl,PredictVect);
                   % [model,accuracy]=fetchOutputs(f(k));
                    trainpredict=model.predictFcn(tbl);
                    trainaccuracy=mean(trainpredict==tbl.currentVEI);

                    Result(iResult).NumberVolcano.morph=height(tbl);
                    Result(iResult).Model.morph=model;
                    Result(iResult).Accuracy.morph=accuracy;
                    Result(iResult).RelativeAccuracy.morph=accuracy-Result(iResult).BaselineAccuracy;
                    Result(iResult).RelativeAccuracyZero.morph=accuracy-Result(iResult).ZeroRule;
                    Result(iResult).TrainAccuracy.morph=trainaccuracy;
                    Result(iResult).TrainRelativeAccuracy.morph=trainaccuracy-Result(iResult).BaselineAccuracy;

            elseif k == 7
                % all numeric called: num
                    tbl=Tbl_VEI_allcat;
                    Imissing=cellfun(@(x) isequal(x,'999'),tbl.RockType);
                    tbl=tbl(~Imissing,:);
                    Imissing=isnan(tbl.eruptiontimediff);
                    tbl=tbl(~Imissing,:);
                    Imissing=isnan(tbl.eruptionduration);
                    tbl=tbl(~Imissing,:);
                    
                    p=tbl.Properties.VariableNames;
                    for i=1:length(p),
                        II(i)=isnumeric(tbl.(i));
                    end
                     I=cellfun(@(x) isequal(x,'currentVEI')|isequal(x,'name'),p);
                    PredictVect=p(~I&II);
                    [model,accuracy] = trainClassifier_Generic(tbl, PredictVect);
                    %f(k)=parfeval(@trainClassifier_Generic,2,tbl,PredictVect);
                    %[model,accuracy]=fetchOutputs(f(k));
                    trainpredict=model.predictFcn(tbl);
                    trainaccuracy=mean(trainpredict==tbl.currentVEI);


                    Result(iResult).NumberVolcano.num=height(tbl);
                    Result(iResult).Model.num=model;
                    Result(iResult).Accuracy.num=accuracy;
                    Result(iResult).RelativeAccuracy.num=accuracy-Result(iResult).BaselineAccuracy;
                    Result(iResult).RelativeAccuracyZero.num=accuracy-Result(iResult).ZeroRule;
                    Result(iResult).TrainAccuracy.num=trainaccuracy;
                    Result(iResult).TrainRelativeAccuracy.num=trainaccuracy-Result(iResult).BaselineAccuracy;

            elseif k == 8
                % Tectonic Setting called: tec
                    tbl=Tbl_VEI_allcat;
                    p=tbl.Properties.VariableNames;
                    I=cellfun(@(x) isequal(x,'TectonicSetting'),p);
                    PredictVect=p(I);
                    [model,accuracy] = trainClassifier_Generic(tbl, PredictVect);
                    %f(k)=parfeval(@trainClassifier_Generic,2,tbl,PredictVect);
                    %[model,accuracy]=fetchOutputs(f(k));
                    trainpredict=model.predictFcn(tbl);
                    trainaccuracy=mean(trainpredict==tbl.currentVEI);

                    Result(iResult).NumberVolcano.tec=height(tbl);
                    Result(iResult).Model.tec=model;
                    Result(iResult).Accuracy.tec=accuracy;
                    Result(iResult).RelativeAccuracy.tec=accuracy-Result(iResult).BaselineAccuracy;
                    Result(iResult).RelativeAccuracyZero.tec=accuracy-Result(iResult).ZeroRule;
                    Result(iResult).TrainAccuracy.tec=trainaccuracy;
                    Result(iResult).TrainRelativeAccuracy.tec=trainaccuracy-Result(iResult).BaselineAccuracy;
            end

        end

        filename='Result';
        save(filename,'Result'); 
        

    end
    end
end

%% Evaluate results for unhistoried volcanoes
% last index of Result is history. Need to use corresponding models for
% each
load Result;

for iVEI=0:3
    %Varying Time thresholds
    for itime = 1:2
            for i_bol = 1:2
                if(~isempty(Result(iVEI+1,itime,i_bol,2).Volcano))
                    Rhistoried=Result(iVEI+1,itime,i_bol,1); % corresponding historied result
                    Runhistoried=Result(iVEI+1,itime,i_bol,2);
                 
% make prediction on unhistoried  & check accuracy

                    tbl1=Runhistoried.AttributeTable.Tbl_VEI_allcat;
                    Imissing=cellfun(@(x) isequal(x,'999'),tbl1.RockType);
                    tbl2=tbl1(~Imissing,:);

                    Predict_tec=Rhistoried.Model.tec.predictFcn(tbl1);
                    Predict_dom=Rhistoried.Model.dom.predictFcn(tbl2);
                    Predict_morph=Rhistoried.Model.morph.predictFcn(tbl1);
                    Predict_cat=Rhistoried.Model.cat.predictFcn(tbl2);
                    
                    Accuracy_tec=mean(Predict_tec'==Runhistoried.currentVEI);
                    Accuracy_dom=mean(Predict_dom'==Runhistoried.currentVEI(~Imissing));
                    Accuracy_morph=mean(Predict_morph'==Runhistoried.currentVEI);
                    Accuracy_cat=mean(Predict_cat'==Runhistoried.currentVEI(~Imissing));

                    Result(iVEI+1,itime,i_bol,2).Accuracy.tec=Accuracy_tec;
                    Result(iVEI+1,itime,i_bol,2).Accuracy.dom=Accuracy_dom;
                    Result(iVEI+1,itime,i_bol,2).Accuracy.morph=Accuracy_morph;
                    Result(iVEI+1,itime,i_bol,2).Accuracy.cat=Accuracy_cat;

                    Result(iVEI+1,itime,i_bol,2).RelativeAccuracy.tec=Accuracy_tec-Result(iVEI+1,itime,i_bol,2).BaselineAccuracy;
                    Result(iVEI+1,itime,i_bol,2).RelativeAccuracy.dom=Accuracy_dom-Result(iVEI+1,itime,i_bol,2).BaselineAccuracy;
                    Result(iVEI+1,itime,i_bol,2).RelativeAccuracy.morph=Accuracy_morph-Result(iVEI+1,itime,i_bol,2).BaselineAccuracy;
                    Result(iVEI+1,itime,i_bol,2).RelativeAccuracy.cat=Accuracy_cat-Result(iVEI+1,itime,i_bol,2).BaselineAccuracy;

                    Result(iVEI+1,itime,i_bol,2).RelativeAccuracyZero.tec=Accuracy_tec-Result(iVEI+1,itime,i_bol,2).ZeroRule;
                    Result(iVEI+1,itime,i_bol,2).RelativeAccuracyZero.dom=Accuracy_dom-Result(iVEI+1,itime,i_bol,2).ZeroRule;
                    Result(iVEI+1,itime,i_bol,2).RelativeAccuracyZero.morph=Accuracy_morph-Result(iVEI+1,itime,i_bol,2).ZeroRule;
                    Result(iVEI+1,itime,i_bol,2).RelativeAccuracyZero.cat=Accuracy_cat-Result(iVEI+1,itime,i_bol,2).ZeroRule;

                    NVolc=length(Runhistoried.currentVEI);
                    Result(iVEI+1,itime,i_bol,2).NumberVolcano.tec=height(tbl1);
                    Result(iVEI+1,itime,i_bol,2).NumberVolcano.dom=height(tbl2);
                    Result(iVEI+1,itime,i_bol,2).NumberVolcano.morph=height(tbl1);
                    Result(iVEI+1,itime,i_bol,2).NumberVolcano.cat=height(tbl2);

                end
            end
    end
end



   filename='Result';
   save(filename,'Result'); 

%% Plot bar graphs 

load Result;

% After 1500
HistoriedAccuracy=[Result(:,1,:,1).Accuracy];
UnHistoriedAccuracy=[Result(:,1,:,2).Accuracy];
HistoriedRelativeAccuracy=[Result(:,1,:,1).RelativeAccuracy]; % can put zero here
UnHistoriedRelativeAccuracy=[Result(:,1,:,2).RelativeAccuracy];
HistoriedNumbers=[Result(:,1,:,1).NumberVolcano];
UnHistoriedNumbers=[Result(:,1,:,2).NumberVolcano];
HistoriedErr=[Result(:,1,:,1).BaselineStd];
UnHistoriedErr=[Result(:,1,:,2).BaselineStd];
    
RelAccYlims=[0 0.33];
RelAccYlimsPos=[0 0.2];
errorbarflag=1;

figure(1)
subplot(3,1,3)
PlotTotalAccuracy;
title('Historied Eruptions After 1500')
legend('off')
text(0.01,0.9,'\bf c','Units','normalized','FontSize',24)

%figure(2)
subplot(3,1,1)
PlotRelativeAccuracy;
title('Historied Eruptions After 1500')
text(0.01,0.9,'\bf a','Units','normalized','FontSize',24)

%figure(3)
%h=bar([[UnHistoriedAccuracy.cat]; [UnHistoriedAccuracy.tec]; [UnHistoriedAccuracy.dom]; [UnHistoriedAccuracy.morph]]);
HN=[[UnHistoriedNumbers.cat];[UnHistoriedNumbers.tec]
    [UnHistoriedNumbers.dom]; [UnHistoriedNumbers.morph]];


subplot(3,1,2)
PlotUnhistRelativeAccuracy;
title('Unhistoried Eruptions After 1500')

% All times

HistoriedAccuracy=[Result(:,2,:,1).Accuracy];
UnHistoriedAccuracy=[Result(:,2,:,2).Accuracy];
HistoriedRelativeAccuracy=[Result(:,2,:,1).RelativeAccuracy];
UnHistoriedRelativeAccuracy=[Result(:,2,:,2).RelativeAccuracy];
HistoriedNumbers=[Result(:,2,:,1).NumberVolcano];
UnHistoriedNumbers=[Result(:,2,:,2).NumberVolcano];

errorbarflag=1;

figure(11)


%figure(12)

subplot(3,1,3)
PlotTotalAccuracy;
legend('off')
title('Historied Eruptions All Times')
text(0.01,0.9,'\bf c','Units','normalized','FontSize',24)

subplot(3,1,1)
PlotRelativeAccuracy;

HN=[[UnHistoriedNumbers.cat];[UnHistoriedNumbers.tec]
     [UnHistoriedNumbers.dom]; [UnHistoriedNumbers.morph]];
title('Historied Eruptions All Times')
text(0.01,0.9,'\bf a','Units','normalized','FontSize',24)


subplot(3,1,2)
PlotUnhistRelativeAccuracy;
title('Unhistoried Eruptions All Times')

% figure(13)
% h=bar([[UnHistoriedAccuracy.cat]; [UnHistoriedAccuracy.tec]; [UnHistoriedAccuracy.dom]; [UnHistoriedAccuracy.morph]]);

% 
% 
% legend([Result(:,2,:,2).VolcanoDatasetName])
% x=[{'All Categories'},{'Tectonic Setting'},{'Dominant Rock Type'},{'Morphology'}];
% set(gca,'XTickLabel',x)
% 
% for ih=1:length(h),
%     text(h(ih).XEndPoints,h(ih).YData,num2cell(HN(:,ih)), ...
%         'HorizontalAlignment','right','VerticalAlignment','middle', ...
%         'Color', 'w', 'Rotation', 90,'Fontsize',11)
% end
% ylabel('Total Accuracy ')
% ylim([0 1])

%%
%Plot Zero-rule relative accuracies
errorbarflag=0;

% After 1500
HistoriedAccuracy=[Result(:,1,:,1).Accuracy];
UnHistoriedAccuracy=[Result(:,1,:,2).Accuracy];
HistoriedRelativeAccuracy=[Result(:,1,:,1).RelativeAccuracyZero]; 
UnHistoriedRelativeAccuracy=[Result(:,1,:,2).RelativeAccuracyZero];
HistoriedNumbers=[Result(:,1,:,1).NumberVolcano];
UnHistoriedNumbers=[Result(:,1,:,2).NumberVolcano];
HistoriedErr=[Result(:,1,:,1).BaselineStd];
UnHistoriedErr=[Result(:,1,:,2).BaselineStd];
    
RelAccYlims=[-0.1 0.1];

figure(31)
PlotTotalAccuracy;
clf;
%figure(2)
subplot(2,1,1)
PlotRelativeAccuracy;
title('Historied Eruptions After 1500')
text(0.01,0.9,'\bf a','Units','normalized','FontSize',24)
ylim([-0.2 0.2])

%figure(3)
%h=bar([[UnHistoriedAccuracy.cat]; [UnHistoriedAccuracy.tec]; [UnHistoriedAccuracy.dom]; [UnHistoriedAccuracy.morph]]);
HN=[[UnHistoriedNumbers.cat];[UnHistoriedNumbers.tec]
    [UnHistoriedNumbers.dom]; [UnHistoriedNumbers.morph]];


subplot(2,1,2)
PlotUnhistRelativeAccuracy;
title('Unhistoried Eruptions After 1500')
ylim([-0.05 0.2])

% All times

HistoriedAccuracy=[Result(:,2,:,1).Accuracy];
UnHistoriedAccuracy=[Result(:,2,:,2).Accuracy];
HistoriedRelativeAccuracy=[Result(:,2,:,1).RelativeAccuracyZero];
UnHistoriedRelativeAccuracy=[Result(:,2,:,2).RelativeAccuracyZero];
HistoriedNumbers=[Result(:,2,:,1).NumberVolcano];
UnHistoriedNumbers=[Result(:,2,:,2).NumberVolcano];


figure(41)
PlotTotalAccuracy;
clf;

subplot(2,1,1)
PlotRelativeAccuracy;

HN=[[UnHistoriedNumbers.cat];[UnHistoriedNumbers.tec]
     [UnHistoriedNumbers.dom]; [UnHistoriedNumbers.morph]];
title('Historied Eruptions All Times')
text(0.01,0.9,'\bf a','Units','normalized','FontSize',24)
ylim([-0.2 0.2])

subplot(2,1,2)
PlotUnhistRelativeAccuracy;
title('Unhistoried Eruptions All Times')
ylim([-0.05 0.2])

%%
% Plot training accuracies

load Result;

% After 1500
HistoriedAccuracy=[Result(:,1,:,1).TrainAccuracy];
UnHistoriedAccuracy=[Result(:,1,:,2).TrainAccuracy];
HistoriedRelativeAccuracy=[Result(:,1,:,1).TrainRelativeAccuracy];
UnHistoriedRelativeAccuracy=[Result(:,1,:,2).TrainRelativeAccuracy];
HistoriedNumbers=[Result(:,1,:,1).NumberVolcano];
UnHistoriedNumbers=[Result(:,1,:,2).NumberVolcano];

RelAccYlims=[-0.15 0.3];
RelAccYlimsPos=[0 0.2];

for j=1:4,
HistoriedAccuracy(j).last=Result(j,1,1,1).Accuracy.last;
HistoriedAccuracy(j).median=Result(j,1,1,1).Accuracy.median;
HistoriedAccuracy(j).mode=Result(j,1,1,1).Accuracy.mode;
HistoriedAccuracy(j).min=Result(j,1,1,1).Accuracy.min;
HistoriedAccuracy(j).max=Result(j,1,1,1).Accuracy.max;

HistoriedRelativeAccuracy(j).last=Result(j,1,1,1).RelativeAccuracy.last;
HistoriedRelativeAccuracy(j).median=Result(j,1,1,1).RelativeAccuracy.median;
HistoriedRelativeAccuracy(j).mode=Result(j,1,1,1).RelativeAccuracy.mode;
HistoriedRelativeAccuracy(j).min=Result(j,1,1,1).RelativeAccuracy.min;
HistoriedRelativeAccuracy(j).max=Result(j,1,1,1).RelativeAccuracy.max;
end
j=5;
HistoriedAccuracy(j).last=Result(3,1,2,1).Accuracy.last;
HistoriedAccuracy(j).median=Result(3,1,2,1).Accuracy.median;
HistoriedAccuracy(j).mode=Result(3,1,2,1).Accuracy.mode;
HistoriedAccuracy(j).min=Result(3,1,2,1).Accuracy.min;
HistoriedAccuracy(j).max=Result(3,1,2,1).Accuracy.max;

HistoriedRelativeAccuracy(j).last=Result(3,1,2,1).RelativeAccuracy.last;
HistoriedRelativeAccuracy(j).median=Result(3,1,2,1).RelativeAccuracy.median;
HistoriedRelativeAccuracy(j).mode=Result(3,1,2,1).RelativeAccuracy.mode;
HistoriedRelativeAccuracy(j).min=Result(3,1,2,1).RelativeAccuracy.min;
HistoriedRelativeAccuracy(j).max=Result(3,1,2,1).RelativeAccuracy.max;

figure(21)

HA=[[HistoriedAccuracy.allcat]; [HistoriedAccuracy.cat]; ...
    [HistoriedAccuracy.dom]; [HistoriedAccuracy.etd]; ...
    [HistoriedAccuracy.ed]; [HistoriedAccuracy.morph]; ...
    [HistoriedAccuracy.num]; [HistoriedAccuracy.tec]; ...
    [HistoriedAccuracy.last]; ...
    [HistoriedAccuracy.median]; ...
    [HistoriedAccuracy.mode]; ...
    [HistoriedAccuracy.min]; ...
    [HistoriedAccuracy.max]];

HN=[[HistoriedNumbers.allcat]; [HistoriedNumbers.cat]; ...
    [HistoriedNumbers.dom]; [HistoriedNumbers.etd]; ...
    [HistoriedNumbers.ed]; [HistoriedNumbers.morph]; ...
    [HistoriedNumbers.num]; [HistoriedNumbers.tec]
    [HistoriedNumbers.Simple]; [HistoriedNumbers.Simple]; ...
    [HistoriedNumbers.Simple]; [HistoriedNumbers.Simple]; ...
    [HistoriedNumbers.Simple]];

h=bar(HA);
legend([Result(:,1,:,1).VolcanoDatasetName])
%x=[{'All Categories'},{'Tectonic Setting'},{'Dominant Rock Type'},{'Morphology'}];
Attributes = {'All Attributes','Categorical Data','Dominant Rock Type', ...
    'Repose Time','Eruption Duration','Morphology','All Numerical Data','Tectonic Setting' ...
    'Last VEI','Median VEI','Mode VEI','Min VEI','Max VEI'};

set(gca,'XTickLabel',Attributes)

for ih=1:length(h),
    text(h(ih).XEndPoints,h(ih).YData,num2cell(HN(:,ih)), ...
        'HorizontalAlignment','right','VerticalAlignment','middle', ...
        'Color', 'w', 'Rotation', 90,'Fontsize',11)
end
%'HorizontalAlignment','right','VerticalAlignment','middle', 'Color', 'w', 'Rotation', 90,'Fontsize',11)
ylabel('Total Accuracy ')
ylim([0 1])
% separate simple predictors
hold on
plot(([1 1]*h(5).XEndPoints(8)+h(1).XEndPoints(9))/2,[0 1],'k--')
hold off

figure(22)
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
legend([Result(:,1,:,1).VolcanoDatasetName])
%x=[{'All Categories'},{'Tectonic Setting'},{'Dominant Rock Type'},{'Morphology'}];
Attributes = {'All Attributes','Categorical Data','Dominant Rock Type', ...
    'Repose Time','Eruption Duration','Morphology','All Numerical Data','Tectonic Setting' ...
    'Last VEI','Median VEI','Mode VEI','Min VEI','Max VEI'};

set(gca,'XTickLabel',Attributes)

for ih=1:length(h),
    text(h(ih).XEndPoints,h(ih).YData,num2cell(HN(:,ih)), ...
        'HorizontalAlignment','right','VerticalAlignment','middle', ...
        'Color', 'w', 'Rotation', 90,'Fontsize',11)
end
% separate simple predictors
hold on
plot(([1 1]*h(5).XEndPoints(8)+h(1).XEndPoints(9))/2,[0 1],'k--')
hold off
ylim(RelAccYlims)


%% Pie chart
% compare distribution of morphologies in full and correctly predicted datasets
%
figure(50)
MakePieChart(Result(3,1,2,1),Result(3,1,2,1).Model.dom,'morphology');
subplot(1,2,1)
text(0.01,0.95,'\bf a','FontSize',18,'Units','normalized')

figure(51)
MakePieChart(Result(3,1,2,1),Result(3,1,2,1).Model.allcat,'currentVEI');
subplot(1,2,1)
text(0.01,0.95,'\bf b','FontSize',18,'Units','normalized')
