%% Eye-movement Pre-analysis
% Input raw data, fixation position, distance threshold
% Return the disqualified trial numbers
% --------------------------------------------
% [InvalidFixTrial,NoSaccadeTrial]=EyelinkPreAnalysis(RawData,FixPos,FixThreshold)
% --------------------------------------------
% RawData is a 5-column array, including:
% Trial index/ Gaze X / Gaze Y/ Saccade index/ Fixation index
%
% MAC lab, ECNU, 2018.11.12

function [InvalidFixTrial,NoSaccadeTrial]=EyelinkPreAnalysis(RawData,FixPos,FixThreshold)
InvalidFixTrial=[];
NoSaccadeTrial=[];
BlockLength=max(RawData(:,1));
for i=1:BlockLength
    CurrentTrialBool= RawData(:,1)==i;
    CurrentTrial=RawData(CurrentTrialBool,:); % Pick out single-trial data
    Bool_1stFix= CurrentTrial(:,5)==1; 
    XY_1stFix=CurrentTrial(Bool_1stFix,:); 
    Aver_1stFixPos=[mean(XY_1stFix(:,2)),mean(XY_1stFix(:,3))];
    % Average position of 1st fixation
    FixDist=sqrt((Aver_1stFixPos(1)-FixPos(1))^2+(Aver_1stFixPos(2)-FixPos(2))^2);
    % Count geometric distance between 1st fixation and the fixation point
    if FixDist>FixThreshold || length(XY_1stFix)<10
        InvalidFixTrial=horzcat(InvalidFixTrial,i);
    end
    if max(CurrentTrial(:,4))<1
        NoSaccadeTrial=horzcat(NoSaccadeTrial,i);
    end
end
% disp('Pre-processing is done.')
