%% Visualization for Behavioral Analysis
% Memory-based attentional capture
% M.A.C lab in ECNU
% 2018.7.5

clf

figure(1);% between-condition RT
box off
errorbar([1,2],[t_RT_match(1) t_RT_dismatch(1)],[std(data.RT_match(:,1)) std(data.RT_dismatch(:,1))],...
    'LineWidth',1.5)
hold on
errorbar([1.1,2.1],[t_RT_match(2) t_RT_dismatch(2)],[std(data.RT_match(:,2)) std(data.RT_dismatch(:,2))],...
    'LineWidth',1.5)

ylabel('Mean RT(s)')
xlabel('Distractor')

set(gca, 'XTick',[1.1,2.1] );
set(gca, 'XTicklabel',{'Match','Dismatch'});
set(gca, 'LineWidth',1.5);
legend('1','2')
figure(2);% Raw SD & mean deviation
subplot(1,2,1)
box off
plot(t_rawSD,'.-','LineWidth',1.5)
ylabel('Raw SD')
xlabel('Memory Load')

set(gca, 'XTick',[1,2] );
set(gca, 'XTicklabel',{'1','2'});
set(gca, 'LineWidth',1.5);
subplot(1,2,2)
box off
plot(t_mDev,'.-','LineWidth',1.5)
ylabel('Mean Deviation')
xlabel('Memory Load')

set(gca, 'XTick',[1,2] );
set(gca, 'XTicklabel',{'1','2'});
set(gca, 'LineWidth',1.5);
if model=='y'
    figure(3);% Model fitting results
    subplot(1,2,1)
    box off
    set(gca, 'LineWidth',1.5);
    plot(t_mSD,'.-','LineWidth',1.5)
    ylabel('Model SD')
    xlabel('Memory Load')
    
    set(gca, 'XTick',[1,2] );
    set(gca, 'XTicklabel',{'1','2'});
    subplot(1,2,2)
    box off
    set(gca, 'LineWidth',1.5);
    plot(t_g,'.-','LineWidth',1.5)
    ylabel('g')
    xlabel('Memory Load')
    
    set(gca, 'XTick',[1,2] );
    set(gca, 'XTicklabel',{'1','2'});
end
figure(4);% ACC
box off
set(gca, 'LineWidth',1.5);
plot(t_acc,'.-','LineWidth',1.5)
ylabel('Memory Test Accuracy')
xlabel('Memory Load')

set(gca, 'XTick',[1,2] );
set(gca, 'XTicklabel',{'1','2'});
figure(5);% Correlation between rawSD/mean deviation/g/kappa & attentional capture
for i=1:2
    subplot(4,3,i)
    box off
    set(gca, 'LineWidth',1.5);
    x1 =data.rawSD(:,i);
    y1 =data.RT_cap(:,i);
    [p,s] = polyfit(x1,y1,1);
    [yfit,dy] = polyconf(p,x1,s,'predopt','curve');
    a=[x1 yfit];
    a=sortrows(a,1);
    x1=a(:,1);
    yfit=a(:,2);
    line(x1,yfit,'color',[0 0 0],'Linewidth',1)
    hold on;
    plot(x1,y1,'.')
end
for i=1:2
    subplot(2,4,i+2)
    set(gca, 'LineWidth',1.5);
    box off
    x1 =data.mDev(:,i);
    y1 =data.RT_cap(:,i);
    [p,s] = polyfit(x1,y1,1);
    [yfit,dy] = polyconf(p,x1,s,'predopt','curve');
    a=[x1 yfit];
    a=sortrows(a,1);
    x1=a(:,1);
    yfit=a(:,2);
    line(x1,yfit,'color',[0 0 0],'Linewidth',1)
    hold on;
    plot(x1,y1,'.')
end
for i=1:2
    subplot(2,4,i+4)
    set(gca, 'LineWidth',1.5);
    box off
    x1 =data.g(:,i);
    y1 =data.RT_cap(:,i);
    [p,s] = polyfit(x1,y1,1);
    [yfit,dy] = polyconf(p,x1,s,'predopt','curve');
    a=[x1 yfit];
    a=sortrows(a,1);
    x1=a(:,1);
    yfit=a(:,2);
    line(x1,yfit,'color',[0 0 0],'Linewidth',1)
    hold on;
    plot(x1,y1,'.')
end
for i=1:2
    subplot(2,4,i+6)
    set(gca, 'LineWidth',1.5);
    box off
    x1 =data.mSD(:,i);
    y1 =data.RT_cap(:,i);
    [p,s] = polyfit(x1,y1,1);
    [yfit,dy] = polyconf(p,x1,s,'predopt','curve');
    a=[x1 yfit];
    a=sortrows(a,1);
    x1=a(:,1);
    yfit=a(:,2);
    line(x1,yfit,'color',[0 0 0],'Linewidth',1)
    hold on;
    plot(x1,y1,'.')
end
if class=='y'
    figure(6); % RT comparisions after classfication
    x_fig6=[1:length(thre)-1; (1:length(thre)-1)+length(thre)-1];
    y_fig6=[t_RT_match_class(:,1)'; t_RT_match_class(:,2)';...
        t_RT_dismatch_class(:,1)'; t_RT_dismatch_class(:,2)'];
    sd_fig6=[sd_RT_match_class(:,1)'; sd_RT_match_class(:,2)';...
        sd_RT_dismatch_class(:,1)'; sd_RT_dismatch_class(:,2)'];
    y_fig6_probe=[t_RT_probematch_class(:,1)'; t_RT_probematch_class(:,2)';...
        t_RT_dismatch_class(:,1)'; t_RT_probedismatch_class(:,2)'];
    sd_fig6_probe=[sd_RT_probematch_class(:,1)'; sd_RT_probematch_class(:,2)';...
        sd_RT_dismatch_class(:,1)'; sd_RT_probedismatch_class(:,2)'];
    box off
    set(gca,'LineWidth',1.5);
    subplot(2,1,1)
    for i=1:2
        hold on
        errorbar(x_fig6(i,:)-0.05,y_fig6(i,:),sd_fig6(i,:)/sqrt(13),'Color',[1,0,0])
        hold on
        errorbar(x_fig6(i,:)+0.05,y_fig6(i+2,:),sd_fig6(i+2,:)/sqrt(13),'Color',[0,0,1])
    end
    hold on
    errorbar(x_fig6(2,:)+5-0.05,y_fig6_probe(2,:),sd_fig6_probe(3,:)/sqrt(13),'Color',[1,0,0])
    errorbar(x_fig6(2,:)+5+0.05,y_fig6(4,:),sd_fig6(4,:)/sqrt(13),'Color',[0,0,1])
    box off
    set(gca,'LineWidth',1.5);
    subplot(2,1,2)
    for i=1:2
        hold on
        plot(x_fig6(i,:),y_fig6_probe(i,:),'Color',[1,0,0])
        hold on
        plot(x_fig6(i,:),y_fig6(i+2,:),'Color',[0,0,1])
    end
    %     subplot(3,1,3)
    %     for i=1:3
    %         hold on
    %         plot(x_fig6(i,:),y_fig6_probe(i+3,:),'Color',[1,0,0])
    %         hold on
    %         plot(x_fig6(i,:),y_fig6(i+3,:),'Color',[0,0,1])
    %     end
    
%     figure(7) % Distribution of probe-matching trials among precision levels
%     box off
%     set(gca,'LineWidth',1.5);
%     mratio=zeros(sortnum, sortnum);
%     mratio_2=mratio;
%     for i=1:sortnum
%         for j=1:SS
%             mratio(i,j)=mean(data.ratio_sort_probematch(:,j,i));
%         end
%     end
%     for i=1:SS
%         mratio_2(i)=mean(data.ratio_probe_match(:,i));
%     end
%     bar(1:3,mratio(1,:)-mratio_2(1),'r')
%     hold on
%     bar([1:3]+4,mratio(2,:)-mratio_2(2),'r')
%     hold on
%     bar([1:3]+8,mratio(3,:)-mratio_2(3),'r')
end

x_fig8=vertcat(x_fig6, 11:15);
%     y_fig8=[t_RT_match_class(:,1)'; t_RT_match_class(:,2)'; t_RT_probematch_class(:,2)';...
%         t_RT_dismatch_class(:,1)'; t_RT_dismatch_class(:,2)'; t_RT_probedismatch_class(:,2)'];
%     sd_fig8=[sd_RT_match_class(:,1)'; sd_RT_match_class(:,2)'; sd_RT_probematch_class(:,2)';...
%         sd_RT_dismatch_class(:,1)'; sd_RT_dismatch_class(:,2)'; sd_RT_probedismatch_class(:,2)'];
y_fig8=vertcat(y_fig6, y_fig6_probe(3:4,:));
sd_fig8=vertcat(sd_fig6, sd_fig6_probe(3:4,:));
for i=1:3
        hold on
        errorbar(x_fig8(i,:)-0.05,y_fig8(i,:),sd_fig8(i,:)/sqrt(13),'Color',[1,0,0])
        hold on
        errorbar(x_fig8(i,:)+0.05,y_fig8(i+3,:),sd_fig8(i+2,:)/sqrt(13),'Color',[0,0,1])
    end
