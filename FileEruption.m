%% Function For Filing Eruption Data 

function eruption=FileEruption(Row)

%Pulls out the EruptionNumber  from the eruption data set and makes it a
%row called number  
eruption.number=Row.EruptionNumber;
%Pulls out the VEI from the eruption data set and makes it a
%row called vei
eruption.vei=Row.VEI;
%Pulls out the StartYear from the eruption data set and makes it a
%row called startyear 
eruption.startyear=Row.StartYear;
%Pulls out the StartMonth from the eruption data set and makes it a
%row called startmonth
eruption.startmonth=Row.StartMonth;
%Pulls out the StartDay from the eruption data set and makes it a
%row called startday
eruption.startday=Row.StartDay;
%Pulls out the EndYear from the eruption data set and makes it a
%row called endyear
eruption.endyear=Row.EndYear;
%Pulls out the EndMonth from the eruption data set and makes it a
%row called endmonth 
eruption.endmonth=Row.EndMonth;
%Pulls out the EndDay from the eruption data set and makes it a
%row called endday
eruption.endday=Row.EndDay;
end


