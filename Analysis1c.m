%% Analysis 1c
% Memory-based attentional capture
% M.A.C lab in ECNU
% 2018.10.8

clear variables

disp('Welcome, fellow.')
disp('Press any key whenever u''re ready')
pause
model=input('Would u like to do the model fitting this time?(y/n)','s');
probe_match=input('Would u like to calculate precision specifically for the mismatching probe this time?(y/n)','s');
class=input('Would u like to classify the data according to the deviation?(y/n)','s');

% Load data
file1=dir('*t1_mix*.mat');
file2=dir('*t2_mix*.mat');
Nsubj=length(file1)/5; % Sample size
SS=2; % Set size
if class=='y'
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
if probe_match=='y'
    data.RT_probe_match=zeros(Nsubj,2);
    data.RT_probe_dismatch=zeros(Nsubj,2);
    data.probe_beta=zeros(Nsubj,1);
    data.probe_ALLbeta=zeros(Nsubj,4);
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
for cond2=1:SS
    search=[];
    memory_answer=[];
    count2=Nsubj;
    % Merge data
    for j=1:Nsubj
        search=[];
        memory_answer=[];
        search_valid=[];
        if cond2==1
            search_m=[];
            for k=0:4
                load(file1(5*j-k).name)
                A_DataMerging1c;
                search=vertcat(search,search_m);
            end
        else
            for k=0:4
                load(file2(5*j-k).name)
                A_DataMerging1c;
                search=vertcat(search,search_m);
            end
        end
        % Skip technical error
        search_valid=search(search(:,9)~=0,:);
        data.search(j,cond2)=mean(search_valid(:,2));
        % Skip erroneous trials
        search_cor=[];
        search_err=[];
        for i=1:length(search_valid)
            if search_valid(i,2)==1
                search_cor=vertcat(search_cor, search_valid(i,:));
            else
                search_err=vertcat(search_err, search_valid(i,:));
            end
        end
        if probe_match=='y' && cond2==2
            search_probe=[];
            for i=1:length(search_cor)
                if search_cor(i,12)~=0
                    search_probe=vertcat(search_probe, search_cor(i,:));
                end
            end
        end
        % Separate conditions
        search_match=[];
        search_dismatch=[];
        for i=1:length(search_cor)
            if search_cor(i,4)==1
                search_match=vertcat(search_match,search_cor(i,:));
            else
                search_dismatch=vertcat(search_dismatch,search_cor(i,:));
            end
        end
        if probe_match=='y'
            search_probe_match=[];
            search_probe_dismatch=[];
            for i=1:length(search_match)
                if search_match(i,12)==1
                    search_probe_match=vertcat(search_probe_match,search_match(i,:));
                else
                    search_probe_dismatch=vertcat(search_probe_dismatch,search_match(i,:));
                end
            end
            mRT_probe_match=mean(search_probe_match(:,3));
            if cond2==1
                mRT_probe_dismatch=NaN;
            else
                mRT_probe_dismatch=mean(search_probe_dismatch(:,3));
            end
            data.ratio_probe_match(j,cond2)=length(search_probe_match(:,3))/length(search_valid);
            data.RT_probe_match(j,cond2)=mRT_probe_match;
            data.RT_probe_dismatch(j,cond2)=mRT_probe_dismatch;
        end
        % Classify the deviation
        if class=='y'
            search_sort=sortrows(search_cor,10);
            n=size(search_sort);
            sortnum=5;
            thre=0;
            for i=1:sortnum
                thre=[thre round(n(1)*i/sortnum)];
            end
            for i=1:sortnum
                sort_range_1=thre(i)+1;
                sort_range_2=thre(i+1);
                search_sort_temp=search_sort(sort_range_1:sort_range_2,:);
                match_ind=find(search_sort_temp(:,4)==1);
                dismatch_ind=setdiff(1:(sort_range_2-sort_range_1+1),match_ind);
                search_sort_match=search_sort_temp(match_ind,:);
                search_sort_dismatch=search_sort_temp(dismatch_ind,:);
                data.RT_match_class(j,cond2,i)=mean(search_sort_match(:,3));
                data.Dev_match_class(j,cond2,i)=mean(search_sort_match(:,10));
                data.Dev_class(j,cond2,i)=mean(search_sort_temp(:,10));
                if model=='y'
                    fit_match=MemFit(search_sort_match(:,9),StandardMixtureModel);
                    data.mSD_match_class(j,cond2,i)=fit_match.posteriorMean(2);
                    data.g_match_class(j,cond2,i)=fit_match.posteriorMean(1);
                    fit_dismatch=MemFit(search_sort_dismatch(:,9),StandardMixtureModel);
                    data.mSD_dismatch_class(j,cond2,i)=fit_dismatch.posteriorMean(2);
                    data.g_dismatch_class(j,cond2,i)=fit_dismatch.posteriorMean(1);
                end
                data.RT_dismatch_class(j,cond2,i)=mean(search_sort_dismatch(:,3));
                data.Dev_dismatch_class(j,cond2,i)=mean(search_sort_dismatch(:,10));
                if probe_match=='y' && cond2==2
                    sort_probematch_ind=find(search_sort_match(:,12)==1);
                    ratio_probematch=length(sort_probematch_ind)/length(search_sort_match(:,12));
                    data.ratio_sort_probematch(j,cond2,i)=ratio_probematch;
                    sort_probedismatch_ind=find(search_sort_match(:,12)~=1);
                    search_sort_probematch=search_sort_match(sort_probematch_ind,:);
                    search_sort_probedismatch=search_sort_match(sort_probedismatch_ind,:);
                    data.RT_match_class_probematch(j,cond2,i)=mean(search_sort_probematch(:,3));
                    data.Dev_match_class_probematch(j,cond2,i)=mean(search_sort_probematch(:,10));
                    data.RT_match_class_probedismatch(j,cond2,i)=mean(search_sort_probedismatch(:,3));
                    data.Dev_match_class_probedismatch(j,cond2,i)=mean(search_sort_probedismatch(:,10));
                end
            end
        end
%          if class=='y'
%             search_sort=sortrows(search_cor,11);
%             n=size(search_sort);
%             thre1=round(n(1)/4);
%             thre2=round(n(1)*2/4);
%             thre3=round(n(1)*3/4);
%             thre=[0 thre1 thre2 thre3 n(1)];
%             for i=1:4
%                 sort_range_1=thre(i)+1;
%                 sort_range_2=thre(i+1);
%                 search_sort_temp=search_sort(sort_range_1:sort_range_2,:);
%                 match_ind=find(search_sort_temp(:,4)==1);
%                 dismatch_ind=setdiff(1:(sort_range_2-sort_range_1+1),match_ind);
%                 search_sort_match=search_sort_temp(match_ind,:);
%                 search_sort_dismatch=search_sort_temp(dismatch_ind,:);
%                 data.RT_match_class(j,cond2,i)=mean(search_sort_match(:,3));
%                 data.RT_dismatch_class(j,cond2,i)=mean(search_sort_dismatch(:,3));
%                 if probe_match=='y'
%                     sort_probematch_ind=find(search_sort_match(:,12)==1);
%                     ratio_probematch=length(sort_probematch_ind)/length(search_sort_match(:,12));
%                     data.ratio_sort_probematch(j,cond2,i)=ratio_probematch;
%                     sort_probedismatch_ind=find(search_sort_match(:,12)~=1);
%                     search_sort_probematch=search_sort_match(sort_probematch_ind,:);
%                     search_sort_probedismatch=search_sort_match(sort_probedismatch_ind,:);
%                     data.RT_match_class_probematch(j,cond2,i)=mean(search_sort_probematch(:,3));
%                     data.RT_match_class_probedismatch(j,cond2,i)=mean(search_sort_probedismatch(:,3));
%                 end
%             end
%         end
        mRT_match=mean(search_match(:,3));
        mRT_dismatch=mean(search_dismatch(:,3));
        mRT_cap=mRT_match-mRT_dismatch;
        mDev=mean(search_valid(:,10));
        [ang_var,circ_sd]=circ_std(pi*search_valid(:,10)/180);
        acc=mean(search_valid(:,11));
        r=corrcoef(search_match(:,3),search_match(:,10));
        if model=='y'
            fit=MemFit(search_valid(:,10),StandardMixtureModel);
        end
        data.RT_match(j,cond2)=mRT_match;
        data.RT_dismatch(j,cond2)=mRT_dismatch;
        data.RT_cap(j,cond2)=mRT_cap;
        data.rawSD(j,cond2)=circ_sd/pi*180;
        data.mDev(j,cond2)=mDev;
        if model=='y'
            data.mSD(j,cond2)=fit.posteriorMean(2);
            data.g(j,cond2)=fit.posteriorMean(1);
        end
        data.ACC(j,cond2)=acc;
        data.cor(j,cond2)=r(2);
        s_indep=size(search_cor);
        INDEP=[search_cor(:,4), search_cor(:,10), search_cor(:,4).*search_cor(:,10), ones(s_indep(1),1)];
        DEP=search_cor(:,3);
        b=regress(DEP,INDEP);
        data.ALLbeta(j,cond2,:)=b;
        data.beta(j,cond2)=b(3);
        count2=count2-1;
        if cond2==2 && probe_match=='y'
            s_indep=size(search_probe);
            INDEP=[search_probe(:,4), search_probe(:,10), search_probe(:,4).*search_probe(:,10), ones(s_indep(1),1)];
            DEP=search_probe(:,3);
            [b,bint,~,~,stats_regress]=regress(DEP,INDEP);
            data.probe_ALLbeta(j,:)=b;
            data.probe_beta(j)=b(3);
        end
        if count2>1
            fprintf('\n%d subjects left in SS %d\n',[count2,cond2])
        else
            fprintf('\n%d subject left in SS %d\n',[count2,cond2])
        end
    end
end
% Overall Statistics
if class=='y'
    for i=1:length(thre)-1
        sd_RT_match_class(i,:)=[std(data.RT_match_class(:,1,i)) std(data.RT_match_class(:,2,i))];
        sd_RT_dismatch_class(i,:)=[std(data.RT_dismatch_class(:,1,i)) std(data.RT_dismatch_class(:,2,i))];
        t_RT_match_class(i,:)=[mean(data.RT_match_class(:,1,i)) mean(data.RT_match_class(:,2,i))];
        t_RT_dismatch_class(i,:)=[mean(data.RT_dismatch_class(:,1,i)) mean(data.RT_dismatch_class(:,2,i))];
    end
    if probe_match=='y'
        for i=1:length(thre)-1
            sd_RT_probematch_class(i,:)=[std(data.RT_match_class_probematch(:,1,i)) std(data.RT_match_class_probematch(:,2,i))];
            t_RT_probematch_class(i,:)=[mean(data.RT_match_class_probematch(:,1,i)) mean(data.RT_match_class_probematch(:,2,i))];
            sd_RT_probedismatch_class(i,:)=[std(data.RT_match_class_probedismatch(:,1,i)) std(data.RT_match_class_probedismatch(:,2,i))];
            t_RT_probedismatch_class(i,:)=[mean(data.RT_match_class_probedismatch(:,1,i)) mean(data.RT_match_class_probedismatch(:,2,i))];
        end
    end
end
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
%t_probe_beta=mean(data.probe_beta);
save(fullfile('Behav_Result3'))
% Visualization
Vis_ans=input('Would u like to visualize the results now?(y/n)','s');
if Vis_ans=='y'
    Visualization_Behav1c
end
% ttests
for i=1:length(thre)-1
    bf.ttest(data.RT_match_class_probematch(:,2,i), data.RT_dismatch_class(:,2,i))
end