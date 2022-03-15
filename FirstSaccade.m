%% FirstSaccade
% Calculate 1st saccade latency
% Capture 1st saccade landing position
% --------------------------------------------
% [Latency_1stSac,HitTrialInd]=FirstSaccade(RawData,AreaOfInterest,Threshold)
% --------------------------------------------
% Let N denote the trial index, then
% RawData is a N*5 matrix, including:
% Trial index/ Gaze X/ Gaze Y/ Saccade index/ Fixation index
% AreaOfInterest is a N*2 matrix, denotes the central position in the area of interest
% Threshold denotes the minimum distance between a qualified 1st saccade landing 
% position and the position of interest, which constructed an area of interest.
% --------------------------------------------
% Latency_1stSac returns a N*1 vector containing 1st saccade latency results
% in number of observations
% HitTrialInd returns a vector representing the index of trials with
% their 1st saccade landing in the area of interest.
% 
% MAC lab, ECNU, 2018.11.12

function [Latency_1stSac,HitTrialInd]=FirstSaccade(RawData,AreaOfInterest,Threshold)
BlockLength=max(RawData(:,1));
Latency_1stSac=zeros(BlockLength,1);
HitTrialInd=[];
for i=1:BlockLength
    CurrentTrialBool= RawData(:,1)==i;
    CurrentTrial=RawData(CurrentTrialBool,:);
    Bool_1stSac= CurrentTrial(:,4)==1;
    CurrentTrialInd=1:length(CurrentTrial);
    Ind_1stSac=CurrentTrialInd(Bool_1stSac);
    Detect_2ndFix=find(CurrentTrial(:,5)==2, 1);
    if isempty(Detect_2ndFix)
        continue
    end
    SacLatency=Ind_1stSac(1);
    Latency_1stSac(i)=SacLatency;
    Ind_LandSac=Ind_1stSac(end)+1;
    LandPos_1stSac=[CurrentTrial(Ind_LandSac,2),CurrentTrial(Ind_LandSac,3)];
    PosOfInterest=AreaOfInterest(i,:);
    SacLandDist=sqrt((LandPos_1stSac(1)-PosOfInterest(1))^2+(LandPos_1stSac(2)-PosOfInterest(2))^2);
    if SacLandDist<Threshold
        HitTrialInd=horzcat(HitTrialInd,i);
    end
end
