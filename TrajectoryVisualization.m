%% Eye Movement Trajectories Visualization
% MAC lab, ECNU, 2018.11.18

close all
clear all
clc

%% Loading data
Eyefile1=dir('*WMbasedAttentionCapture1*.txt');
Eyefile2=dir('*WMbasedAttentionCapture2*.txt');
Behavfile1=dir('*t1_mix*.mat');
Behavfile2=dir('*t2_mix*.mat');
Nsubj=length(Eyefile1)/5;
for SS=1:2
    tc=1;
    for j=1:Nsubj
        search_m=[];
        memory_answer=[];
        search_valid=[];
        BehavInd=zeros(64,5);
        probe_match='n';
        if SS==1
            search=[];
            for k=0:4
                load(Behavfile1(5*j-k).name)
                A_DataMerging1c_EyeMovement;
                search=vertcat(search,search_m);
            end
        else
            search=[];
            for k=0:4
                load(Behavfile2(5*j-k).name)
                A_DataMerging1c_EyeMovement;
                search=vertcat(search,search_m);
            end
        end
        SearchErrInd_t= find(search(:,2)==0);
        Bool_SearchCor= search(:,2)==1;
        search_cor=search(Bool_SearchCor,:);
        search_err=search(SearchErrInd_t,:);
        search_sort=sortrows(search_cor,10);
        search_sort0=vertcat(search_sort,search_err);
        n=size(search_sort);
        sortnum=3;
        thre=0;
        for i=1:sortnum
            thre=[thre round(n(1)*i/sortnum)];
        end
        Behav_Tag=[];
        for i=1:sortnum
            Behav_Tag=[Behav_Tag; i*ones(thre(i+1)-thre(i),1)];
        end
        Behav_Tag=[Behav_Tag; zeros(length(SearchErrInd_t),1)];
        search_sort2=[search_sort0 Behav_Tag];
        search_sort3=sortrows(search_sort2,12);
        for i=1:5
            BehavInd(:,i)=search_sort3(64*(i-1)+1:64*i,13);
        end
        for k=1:5
            %%
            if SS==1
                EyeData=EyelinkSuperReadTXT(Eyefile1(tc).name);
                load(Behavfile1(tc).name);
                tc=tc+1;
            else
                EyeData=EyelinkSuperReadTXT(Eyefile2(tc).name);
                load(Behavfile2(tc).name)
                tc=tc+1;
            end
            
            %% Omit the empty blocks
            if isempty(EyeData)
                continue
            end
            %% Converting behavior data
            BehavInd=zeros(64,5);
            probe_match='n';
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
            SearchErrInd=find(BehavInd(:,k)==0);
            %% Eye Movement Pre-analysis
            [ErrFix,NoSac]=EyelinkPreAnalysis(EyeData,FixPos,DistThreshold);
            ErrTrialInd=union(ErrFix,NoSac);
            ErrTrialInd=union(ErrTrialInd,SearchErrInd);
            Bool_match= CondInd==1;
            Bool_dismatch= CondInd==0;
            [Bool_match(ErrTrialInd),Bool_dismatch(ErrTrialInd)]=deal(0);
            Ind_match=find(Bool_match==1);
            Ind_dismatch=find(Bool_dismatch==1);
            %% Visualize Trajectories
            dis_ind_match=dis_ind(Ind_match);
            for i=1:length(Ind_match)
                Traj=EyeData(EyeData(:,1)==Ind_match(i),2:3);
                RotAng=pi/4*(dis_ind_match(i)-1)+pi/8;
                Traj2=CRT_2D(Traj,[512,384],RotAng);
                plot(Traj2(:,1),Traj2(:,2),'.','Color',[192,192,192])
                hold on
            end
        end
        fprintf('%dsubjects left in SS%d.\r',Nsubj-j,SS)
    end
end

