%% EyelinkSuperReadTXT
% MAC lab, ECNU, 2018.11.13

function [RawData]=EyelinkSuperReadTXT(TXT, Column)
DataID=fopen(TXT);
Data=textscan(DataID,'Trial:%f %f %s %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',2,'TreatAsEmpty',{'NaN'},'EmptyValue',0);
fclose(DataID);
RawData=zeros(length(Data{Column(1)}), length(Column));
for i=1:length(Column)
    RawData(:,i)=Data{Column(i)};
end
end

