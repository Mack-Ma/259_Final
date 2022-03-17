%% Analysis 1c
% Memory-based attentional capture
% M.A.C lab in ECNU
% 2018.10.8
% Revised by Tianye Ma, 3/14/2022

clear variables
close all
clc

disp('Hahahahaha, Welcome.')
disp('Press any key whenever u''re ready')
pause
model=input('Would u like to fit the models this time?(y/n)','s');
classify=input('Would u like to classify the data according to the deviations?(y/n)','s');

% Test required toolboxes
if exist("MemFit",'file')==0 && isequal(model, 'y')
    error('Memtoolbox is needed for the analysis')
end
if exist("circ_std",'file')==0
    error('CircStats toolbox is needed for the analysis')
end

% Load data
addpath('Data_Behav'); % add datafile path
subj_list=[1 2 3 4 6]; % A list of valid subjects
SS_list=[1 2]; % set size
block_list=1:2:9; % block indices
SS=length(SS_list); % N Set size
Nsubj=length(subj_list); % N subjects 

% Pre-allocating space
if isequal(classify,'y')
    Nclass=3;
    data.RT_match_class=zeros(Nsubj,2,Nclass);
    data.RT_dismatch_class=zeros(Nsubj,2,Nclass);
    data.Dev_match_class=zeros(Nsubj,2,Nclass);
    data.Dev_dismatch_class=zeros(Nsubj,2,Nclass);
    data.mSD_match_class=zeros(Nsubj,2,Nclass);
    data.g_match_class=zeros(Nsubj,2,Nclass);
    data.mSD_dismatch_class=zeros(Nsubj,2,Nclass);
    data.g_dismatch_class=zeros(Nsubj,2,Nclass);
    if probe_match=='y'
        data.RT_match_class_probematch=zeros(Nsubj,2,Nclass);
        data.RT_match_class_probedismatch=zeros(Nsubj,2,Nclass);
    end
end
data.search=zeros(Nsubj,2);
data.RT_match=zeros(Nsubj,2);
data.RT_dismatch=zeros(Nsubj,2);
data.RT_cap=zeros(Nsubj,2);
data.rawSD=zeros(Nsubj,2);
data.mDev=zeros(Nsubj,2);
data.mSD=zeros(Nsubj,2);
data.g=zeros(Nsubj,2);
data.ACC=zeros(Nsubj,2);
data.cor=zeros(Nsubj,2);
data.ALLbeta=zeros(Nsubj,2,4);
data.beta=zeros(Nsubj,2);

% Now start analysis
for cond2=1:SS
    count2=Nsubj;
    % Merge data
    for j=1:Nsubj
        search=[];
        memory_answer=[];
        for k=block_list % merge blocks
            filename=['E1c_sub',num2str(subj_list(j)),'_',num2str(k+cond2-1),'_t',num2str(cond2),'_mix.mat'];
            load(filename)
            A_DataMerging1c_EyeMovement;
            search=vertcat(search,search_m);
        end
        % Skip technical error
        search_valid=search(search.mem_ans(:,2)~=0,:); % there is a memory report
        data.search(j,cond2)=mean(search_valid.search_result(:,2)); % search accuracy
        % Skip erroneous trials
        search_cor=search_valid(search_valid.search_result(:,2)==1,:);
        search_err=search_valid(search_valid.search_result(:,2)==0,:);
        % Separate conditions
        search_match=search_cor(search_cor.seq_con==1,:);
        search_dismatch=search_cor(search_cor.seq_con==0,:);
        
        if model=='y'
            save_model=input('Would u like to save the preprocessed dataset before model fitting?(y/n)','s');
            if isequal(save_model,'y')
                save('preprocessed','search_cor')
            end
        end
        
        % Classify the deviation
        if classify=='y'
            search_sort=sortrows(search_cor,"memory_dev");
            n=size(search_sort);
            sortnum=5; % Number of clusters
            thre=zeros(sortnum+1,1);
            for i=2:sortnum+1
                thre(i)=round(n(1)*i/sortnum);
            end
            
            if model=='y'
                save_model1=input('Would u like to save the preprocessed dataset before model fitting?(y/n)','s');
                if isequal(save_model1,'y')
                    save('preprocessed','search_sort')
                end
            end
            
            for i=1:sortnum
                sort_range_1=thre(i)+1;
                sort_range_2=thre(i+1);
                search_sort_temp=search_sort(sort_range_1:sort_range_2,:); % cluster the data
                match_ind=find(search_sort_temp(:,4)==1);
                dismatch_ind=setdiff(1:(sort_range_2-sort_range_1+1),match_ind);
                search_sort_match=search_sort_temp(match_ind,:);
                search_sort_dismatch=search_sort_temp(dismatch_ind,:);
                data.RT_match_class(j,cond2,i)=mean(search_sort_match.search_result(:,3));
                data.Dev_match_class(j,cond2,i)=mean(search_sort_match.memory_dev);
                data.Dev_class(j,cond2,i)=mean(search_sort_temp.memory_dev);
                
                if model=='y' % Mixture modeling?
                    fit_match=MemFit(search_sort_match.memory_dev,StandardMixtureModel);
                    data.mSD_match_class(j,cond2,i)=fit_match.posteriorMean(2);
                    data.g_match_class(j,cond2,i)=fit_match.posteriorMean(1);
                    fit_dismatch=MemFit(search_sort_dismatch.memory_dev,StandardMixtureModel);
                    data.mSD_dismatch_class(j,cond2,i)=fit_dismatch.posteriorMean(2);
                    data.g_dismatch_class(j,cond2,i)=fit_dismatch.posteriorMean(1);
                end
                
                data.RT_dismatch_class(j,cond2,i)=mean(search_sort_dismatch.search_result(:,3));
                data.Dev_dismatch_class(j,cond2,i)=mean(search_sort_dismatch.memory_dev);
            end
        end
        
        [ang_var,circ_sd]=circ_std(pi*search_valid.memory_dev/180); % get circular SDs
        r=corrcoef(search_match.search_result(:,3),search_match.memory_dev);
        data.RT_match(j,cond2)=mean(search_match.search_result(:,3));
        data.RT_dismatch(j,cond2)=mean(search_dismatch.search_result(:,3));
        data.RT_cap(j,cond2)=mean(search_match.search_result(:,3))-...
            mean(search_dismatch.search_result(:,3)); % mean RT difference/capture strength
        data.rawSD(j,cond2)=circ_sd/pi*180; % raw circular SD
        data.mDev(j,cond2)=mean(search_valid.memory_dev); % mean memory deviation
        data.ACC(j,cond2)=mean(search_valid.test_ans); % acc for memory tests;
        data.cor(j,cond2)=r(2);
        s_indep=size(search_cor);
        % GLMs
        INDEP=[search_cor.seq_con, search_cor.memory_dev, search_cor.seq_con.*search_cor.memory_dev, ones(s_indep(1),1)];
        DEP=search_cor.search_result(:,3); % RT
        b=regress(DEP,INDEP);
        data.ALLbeta(j,cond2,:)=b;
        data.beta(j,cond2)=b(3);
        count2=count2-1;
        if model=='y'
            save_model2=input('Would u like to save the preprocessed dataset before model fitting?(y/n)','s');
            if isequal(save_model2,'y')
                save('analysis_modelfree')
            end
            fit=MemFit(search_valid(:,10),StandardMixtureModel);
            data.mSD(j,cond2)=fit.posteriorMean(2);
            data.g(j,cond2)=fit.posteriorMean(1);
        end
        if count2>1
            fprintf('\n%d subjects left in SS %d\n',[count2,cond2])
        else
            fprintf('\n%d subject left in SS %d\n',[count2,cond2])
        end
    end
end
% Get group-level means
t_search=mean(data.search,1);
t_r=[mean(data.cor(:,1)) mean(data.cor(:,2))];
t_RT_match=[mean(data.RT_match(:,1)) mean(data.RT_match(:,2))];
t_RT_dismatch=[mean(data.RT_dismatch(:,1)) mean(data.RT_dismatch(:,2))];
t_mSD=[mean(data.mSD(:,1)) mean(data.mSD(:,2))];
t_g=[mean(data.g(:,1)) mean(data.g(:,2))];
t_acc=[mean(data.ACC(:,1)) mean(data.ACC(:,2))];
t_rawSD=[mean(data.rawSD(:,1)) mean(data.rawSD(:,2))];
t_mDev=[mean(data.mDev(:,1)) mean(data.mDev(:,2))];
t_beta=[mean(data.beta(:,1)), mean(data.beta(:,1))];

save_result=input('Do we save the results?(y/n)','s');
if isequal(save_result,'y')
    save(fullfile('Behav_Result_Hahaha'))
end
% Visualization
Vis_ans=input('Would u like to store the results for visualization now?(y/n)','s');
if isequal(Vis_ans,'y')
    data_Vis=rmfield(data,'ALLbeta');
    writetable(struct2table(data_Vis),'Behav.csv')
end
