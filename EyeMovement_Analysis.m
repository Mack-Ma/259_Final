%% Eye-Movement Analysis_WM-based Attentional Capture
% MAC lab, ECNU, 2018.11.13

close all
clear variables
clc

%% Loading data
Eyefile1=dir('*WM_Capture1*.txt');
Eyefile2=dir('*WM_Capture2*.txt');
Behavfile1=dir('*t1_mix*.mat');
Behavfile2=dir('*t2_mix*.mat');
Nsubj=length(Eyefile1)/5;
EyeQualifiedTrial=zeros(Nsubj,2,2);
TargetTrial=zeros(Nsubj,2,2);
DistrTrial=zeros(Nsubj,2,2);
NontarTrial=zeros(Nsubj,2,2);
Latency=zeros(Nsubj,2,2);
EyeBeta=zeros(Nsubj,2);
EyeAllBeta=zeros(Nsubj,2,4);
COI=[2, 4, 5, 7, 6]; % TrialInd, GazeX, GazeY, Saccade, Fixation
AllTrial_Bool=ones(64,1);
for SS=1:2
    tc=1;
    for j=1:Nsubj
        SacLanding=zeros(3,2);
        TrialNum=zeros(1,2);
        Lat_match=[]; Lat_dismatch=[];
        SacLatency=zeros(1,2);
        Mat=zeros(64*5,5);
        for k=1:5
            if SS==1
                EyeData=EyelinkSuperReadTXT(Eyefile1(tc).name, COI);
                load(Behavfile1(tc).name);
                tc=tc+1;
            else
                EyeData=EyelinkSuperReadTXT(Eyefile2(tc).name, COI);
                load(Behavfile2(tc).name)
                tc=tc+1;
            end
            %% Omit the empty blocks
            if isempty(EyeData)
                continue
            end
            %% Converting behavior data
            FixPos=fix_coord;
            PosArray=search_pos(:,1:2)+r_item;
            TarPos=PosArray(seq_tar_pos,:);
            count=zeros(1,8);
            dis_ind=zeros(64,1);
            for trial=1:64
                if     seq_tar_pos(trial)==1
                    count(1)=count(1)+1;
                    dis_ind(trial)=seq_dis_pos(count(1))+4;
                elseif seq_tar_pos(trial)==2
                    count(2)=count(2)+1;
                    dis_ind(trial)=seq_dis_pos(count(2))+4;
                elseif seq_tar_pos(trial)==3
                    count(3)=count(3)+1;
                    dis_ind(trial)=seq_dis_pos(count(3))+4;
                elseif seq_tar_pos(trial)==4
                    count(4)=count(4)+1;
                    dis_ind(trial)=seq_dis_pos(count(4))+4;
                elseif seq_tar_pos(trial)==5
                    count(5)=count(5)+1;
                    dis_ind(trial)=seq_dis_pos(count(5));
                elseif seq_tar_pos(trial)==6
                    count(6)=count(6)+1;
                    dis_ind(trial)=seq_dis_pos(count(6));
                elseif seq_tar_pos(trial)==7
                    count(7)=count(7)+1;
                    dis_ind(trial)=seq_dis_pos(count(7));
                elseif seq_tar_pos(trial)==8
                    count(8)=count(8)+1;
                    dis_ind(trial)=seq_dis_pos(count(8));
                end
            end
            DisPos=PosArray(dis_ind,:);
            OtherPos=zeros(64,12);
            for trial=1:64
                OtherInd= 1:8~=dis_ind(trial) & 1:8~=seq_tar_pos(trial);
                t_OtherPos=reshape(PosArray(OtherInd,:)',1,12);
                OtherPos(trial,:)=t_OtherPos;
            end
            CondInd=seq_con;
            DistThreshold=sqrt((fix_coord(1)-DisPos(1))^2+(fix_coord(2)-DisPos(2))^2);
            SearchErr_Bool=search_result(:,2);
            SearchErr=find(search_result(:,2)~=1);
            WeirdErr=find(ref_answer(:,1)==0);
            BehavErrTrialInd=union(SearchErr,WeirdErr);
            MemoryDev=2*abs(sub_answer(:,1)-ref_answer(:,1));
            MemoryDev(MemoryDev>180)=360-MemoryDev(MemoryDev>180);
            
            %% 1st Saccade
            [ErrFix,NoSac]=EyelinkPreAnalysis(EyeData,FixPos,DistThreshold);
            EyeErrTrialInd=union(ErrFix,NoSac);
            ErrTrialInd=union(EyeErrTrialInd,BehavErrTrialInd);
            ValidTrialNum=64-length(ErrTrialInd);
            ValidTrialInd=setdiff(1:64,ErrTrialInd);
            ValidTrial_Bool=AllTrial_Bool;
            ValidTrial_Bool(ErrTrialInd)=0;
            [Lat,TarTrialInd]=FirstSaccade(EyeData,TarPos,84);
            [~,DisTrialInd]=FirstSaccade(EyeData,DisPos,84);
            OtherTrialInd0=zeros(64,6);
            for t=1:6
                [~,t_OtherTrialInd]=FirstSaccade(EyeData,OtherPos(:,(2*t-1):2*t),84);
                OtherTrialInd0(t_OtherTrialInd,t)=1;
            end
            OtherTrialInd=[];
            for t=1:6
                OtherTrialInd=union(OtherTrialInd,find(OtherTrialInd0(:,t)==1));
            end
            Bool_match= CondInd==1;
            Bool_dismatch= CondInd==0;
            [Bool_match(ErrTrialInd),Bool_dismatch(ErrTrialInd)]=deal(0);
            Ind_match=find(Bool_match==1);
            Ind_dismatch=find(Bool_dismatch==1);
            Ntrial_match=length(Ind_match);
            Ntrial_dismatch=length(Ind_dismatch);
            Lat_match=vertcat(Lat_match,Lat(Bool_match));
            Lat_dismatch=vertcat(Lat_dismatch,Lat(Bool_dismatch));
            TarTrialInd=setdiff(TarTrialInd,ErrTrialInd);
            TarTrialInd_match=intersect(TarTrialInd,Ind_match);
            TarTrialInd_dismatch=intersect(TarTrialInd,Ind_dismatch);
            DisTrialInd=setdiff(DisTrialInd,ErrTrialInd);
            DisTrialInd_match=intersect(DisTrialInd,Ind_match);
            DisTrialInd_dismatch=intersect(DisTrialInd,Ind_dismatch);
            OtherTrialInd=setdiff(OtherTrialInd,ErrTrialInd);
            OtherTrialInd_match=intersect(OtherTrialInd,Ind_match);
            OtherTrialInd_dismatch=intersect(OtherTrialInd,Ind_dismatch);
            
            TrialNum(1)=TrialNum(1)+Ntrial_match;
            TrialNum(2)=TrialNum(2)+Ntrial_dismatch;
            SacLanding(1,1)=SacLanding(1,1)+length(TarTrialInd_match);
            SacLanding(1,2)=SacLanding(1,2)+length(TarTrialInd_dismatch);
            SacLanding(2,1)=SacLanding(2,1)+length(DisTrialInd_match);
            SacLanding(2,2)=SacLanding(2,2)+length(DisTrialInd_dismatch);
            SacLanding(3,1)=SacLanding(3,1)+length(OtherTrialInd_match);
            SacLanding(3,2)=SacLanding(3,2)+length(OtherTrialInd_dismatch);
            
            Mat((64*(k-1)+1):64*k,:)=[(seq_con+1)', MemoryDev, (seq_con+1)'.*MemoryDev, Lat, ValidTrial_Bool];
            
        end
        Mat=Mat(Mat(:,5)==1,:);
        INDEP=Mat(:,[1, 2, 3, 5]);
        DEP=Mat(:,4);
        b=regress(DEP,INDEP);
        EyeAllBeta(j,SS,:)=b;
        EyeBeta(j,SS)=b(3);
        TargetTrial(j,1,SS)=SacLanding(1,1);
        TargetTrial(j,2,SS)=SacLanding(1,2);
        DistrTrial(j,1,SS)=SacLanding(2,1);
        DistrTrial(j,2,SS)=SacLanding(2,2);
        NontarTrial(j,1,SS)=SacLanding(3,1);
        NontarTrial(j,2,SS)=SacLanding(3,2);
        EyeQualifiedTrial(j,1,SS)=TrialNum(1);
        EyeQualifiedTrial(j,2,SS)=TrialNum(2);
        Latency(j,1,SS)=mean(Lat_match);
        Latency(j,2,SS)=mean(Lat_dismatch);
        fprintf('%dsubjects left in SS%d.\r',Nsubj-j,SS)
    end
end
Psac_Tar=TargetTrial./EyeQualifiedTrial;
Psac_Dis=DistrTrial./EyeQualifiedTrial;
Psac_Nontar=NontarTrial./EyeQualifiedTrial;
%%
save(fullfile(['Behav_Result3' date]))
