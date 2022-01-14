function [Tbl_VEI_allcat] = Create_table_attributes_ED(Volcano,HistoryFlag)

%This code produces all of the Tables for each VEI and time threshold

  
    %Find the unique rock types names
%     [rocktype,~,rocktype_indx] = unique([Volcano.rocktype]','stable');
%     RockType = [rocktype]';

    %Assigning values to the rock types inorder to make them numerical values 
    %This is based on Silica content. 
    %rocktype_values = [57.5,46.5,73,63.3,70,NaN,45,56.2,55,44.75,NaN,49];

%     %Assign the Silica  value to the
%     %various rock types 
%     for r_num = 1:length(rocktype_values)
%        rock_1 = find(rocktype_indx == r_num);
%        rocktype_indx(rock_1,:) = rocktype_values(r_num);
%     end

    % Make change to the Array of the Volcano Morphologies to consoladate the
    % morphologies

    %Finding the Unqiue Morphology of Volcanoes
    [primaryvolcanotype,~,primaryvolcanotype_indx] = unique([Volcano.primaryvolcanotype]','stable');
    PrimaryVolcanoType = categorical(primaryvolcanotype);
    
    % NaN out unknown rocktypes
    for q = 1:length(Volcano)
        if (Volcano(q).rocktype == "")|(Volcano(q).rocktype=="No Data (checked)")
                     Volcano(q).rocktype='999';
        end
    end

    
    for q = 1:length(Volcano)
        if Volcano(q).primaryvolcanotype  == "Complex(es)"
            Volcano(q).primaryvolcanotype = 'Stratovolcano';

        elseif Volcano(q).primaryvolcanotype  == "Complex"
            Volcano(q).primaryvolcanotype = 'Stratovolcano';

        elseif Volcano(q).primaryvolcanotype  == "Compound"
            Volcano(q).primaryvolcanotype = 'Stratovolcano';

        elseif Volcano(q).primaryvolcanotype  == "Cone(s)"
            Volcano(q).primaryvolcanotype = 'Field';

        elseif Volcano(q).primaryvolcanotype  == "Stratovolcano(es)"
            Volcano(q).primaryvolcanotype = 'Stratovolcano';

        elseif Volcano(q).primaryvolcanotype  == "Crater rows"
            Volcano(q).primaryvolcanotype = 'Field';

        elseif Volcano(q).primaryvolcanotype  == "Explosion crater(s)"
            Volcano(q).primaryvolcanotype = 'Field';

        elseif  Volcano(q).primaryvolcanotype  == "Fissure vent"
            Volcano(q).primaryvolcanotype = 'Shield';

        elseif Volcano(q).primaryvolcanotype  == "Fissure vent(s)"
            Volcano(q).primaryvolcanotype = 'Shield';

        elseif Volcano(q).primaryvolcanotype  == "Lava cone"
            Volcano(q).primaryvolcanotype = 'Stratovolcano';

        elseif Volcano(q).primaryvolcanotype  == "Lava dome"
            Volcano(q).primaryvolcanotype = 'Stratovolcano';

        elseif Volcano(q).primaryvolcanotype  == "Lava dome(s)"
            Volcano(q).primaryvolcanotype = 'Stratovolcano';

        elseif Volcano(q).primaryvolcanotype  == "Maar"
            Volcano(q).primaryvolcanotype = 'Field';

        elseif Volcano(q).primaryvolcanotype  == "Maar(s)"
            Volcano(q).primaryvolcanotype = 'Field';

        elseif Volcano(q).primaryvolcanotype  == "Pyroclastic cone"
            Volcano(q).primaryvolcanotype = 'Field';

        elseif Volcano(q).primaryvolcanotype  == "Pyroclastic cone(s)"
            Volcano(q).primaryvolcanotype = 'Field';

        elseif Volcano(q).primaryvolcanotype  == "Pyroclastic shield"
            Volcano(q).primaryvolcanotype = 'Caldera';

        elseif Volcano(q).primaryvolcanotype  == "Caldera(s)"
            Volcano(q).primaryvolcanotype = 'Caldera';

        elseif Volcano(q).primaryvolcanotype  == "Shield(s)"
            Volcano(q).primaryvolcanotype = 'Shield';

        elseif Volcano(q).primaryvolcanotype  == "Stratovolcano?"
            Volcano(q).primaryvolcanotype = 'Stratovolcano';

        elseif Volcano(q).primaryvolcanotype  == "Tuff cone(s)"
            Volcano(q).primaryvolcanotype = 'Field';

        elseif Volcano(q).primaryvolcanotype  == "Tuff ring(s)"
            Volcano(q).primaryvolcanotype = 'Field';

        elseif Volcano(q).primaryvolcanotype  == "Volcanic field"
            Volcano(q).primaryvolcanotype = 'Field';

        end
    end
    if (HistoryFlag)
    for k=1:length(Volcano)

        if (length(Volcano(k).eruption)>1)

   
            %Pulling out the values 
            
            Volcano(k).endyear1 = [Volcano(k).eruption(1).endyear];
            Volcano(k).endmonth1 = [Volcano(k).eruption(1).endmonth];
            Volcano(k).endday1 = [Volcano(k).eruption(1).endday];
            
            Volcano(k).startyear2 = [Volcano(k).eruption(2).startyear];
            Volcano(k).startmonth2 = [Volcano(k).eruption(2).startmonth];
            Volcano(k).startday2 = [Volcano(k).eruption(2).startday];
            
            Volcano(k).endyear2 = [Volcano(k).eruption(2).endyear];
            Volcano(k).endmonth2 = [Volcano(k).eruption(2).endmonth];
            Volcano(k).endday2 = [Volcano(k).eruption(2).endday];
            
            %Calculating the difference and Convert all Categories to years 
            %Eruption time elapsed between eruptions 
            Volcano(k).time2 = datenum(Volcano(k).startyear2,Volcano(k).startmonth2,Volcano(k).startday2);
            Volcano(k).time3 = datenum(Volcano(k).endyear1,Volcano(k).endmonth1,Volcano(k).endday1);
            Volcano(k).time4 = datenum(Volcano(k).endyear2,Volcano(k).endmonth2,Volcano(k).endday2);
            
 
            %Eruption Duration 
            Volcano(k).eruptionduration = Volcano(k).time4 - Volcano(k).time2;
            %Eruption time difference 
            Volcano(k).eruptiontimediff = Volcano(k).time3 - Volcano(k).time2;
 
            %The Silica content of the Rock types
            %Volcano(k).silica = rocktype_indx(k);
 
 
        else
     
 
            %Calculating Eruption dynamics
            Volcano(k).endyear1 = NaN;
            Volcano(k).startyear2 = NaN;
            Volcano(k).endmonth1 = NaN;
            Volcano(k).startmonth2 = NaN;
            Volcano(k).endday1 = NaN;
            Volcano(k).startday2 = NaN;
 
            Volcano(k).startyear1 = NaN;
            Volcano(k).startmonth1 = NaN;
            Volcano(k).startday1 = NaN;
 
            %Datenum varibles 
            Volcano(k).time1 = NaN;
            Volcano(k).time2 = NaN;
            Volcano(k).time3 = NaN;
 
            %Eruption Dynamic Attributes 
            Volcano(k).eruptiontimediff=NaN;
            Volcano(k).eruptionduration=NaN;
 
            %Rocktype as numerical Variable 
            %Volcano(k).silica=NaN;
 
            %The test Variable 
            Volcano(k).currentVEI=NaN;
 
        end
    end
 
    %Removes the [NaN,NaN,NaN] from the datenum varible so the table
    %variable can be concatenated 
    for k = 1:length(Volcano)
        if sum(isnan(Volcano(k).eruptiontimediff))  == 3
            ['WARNING']
            Volcano(k).eruptiontimediff = NaN;
 
        end
 
        if sum(isnan(Volcano(k).eruptionduration))  == 3
 
            Volcano(k).eruptionduration = NaN;
 
        end
 
        %Removes one eruption duration error in Data Base 
        if Volcano(k).eruptionduration < 0 
 
            Volcano(k).eruptionduration = NaN;
 
        end    
    end


    time_thresh = ["after_1500","all_time"];

    % Creating Table of Attributes for Machine Learning Application
 
            I=find(isfinite([Volcano.currentVEI]));
            Volcano=Volcano(I);
            for i=1:length(Volcano),
                Volcano(i).medianVEI=Volcano(i).SimplePredictors.median;
                Volcano(i).maxVEI=Volcano(i).SimplePredictors.max;
                Volcano(i).lastVEI=Volcano(i).SimplePredictors.last;
                Volcano(i).modeVEI=Volcano(i).SimplePredictors.mode;
                Volcano(i).minVEI=Volcano(i).SimplePredictors.min;
            end
            Tbl_VEI_allcat=table([Volcano.currentVEI]',[Volcano.medianVEI]', ...
                [Volcano.primaryvolcanotype]',[Volcano.rocktype]', ...
                [Volcano.tectonicsetting]',[Volcano.lastVEI]', ...
                [Volcano.minVEI]',[Volcano.maxVEI]',[Volcano.modeVEI]', ...
                [Volcano.eruptiontimediff]',[Volcano.eruptionduration]', ...
                [Volcano.name]', ...
                'VariableNames',{'currentVEI','medianVEI', ...
                'morphology', 'RockType', ...
                'TectonicSetting','lastVEI', ...
                'minVEI','maxVEI', 'modeVEI', ...
                'eruptiontimediff','eruptionduration', ...
                'name'});
    else
        Tbl_VEI_allcat=table([Volcano.currentVEI]', ...
                [Volcano.primaryvolcanotype]',[Volcano.rocktype]', ...
                [Volcano.tectonicsetting]', ...
                [Volcano.name]', ...
                'VariableNames',{'currentVEI', ...
                'morphology', 'RockType', ...
                'TectonicSetting', ...
                'name'});
    end
    

