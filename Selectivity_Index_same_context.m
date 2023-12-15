
clc;

clearvars -except stat_p_all resp_all_b resp_all_s resp_all_d resp_all_r

pval = 0.01;
sz = 1;
topcri = 5; % 10 for 10%, 5 for 20%
segno = 5;

edges1 = [-1.05:0.1:1.05];
edges2 = [-1:0.1:1];

s_sig = find(stat_p_all(:,1)<pval);
d_sig = find(stat_p_all(:,3)<pval);
r_sig = find(stat_p_all(:,5)<pval);
s_bsig = find(stat_p_all(:,1)<pval & stat_p_all(:,2)<pval);
s_esig = find(stat_p_all(:,1)<pval & stat_p_all(:,2)>=pval);
s_csig = find(stat_p_all(:,1)<pval & stat_p_all(:,2)>=pval & (stat_p_all(:,4)<pval | stat_p_all(:,6)<pval));
s_nsig = find(stat_p_all(:,1)<pval & stat_p_all(:,2)>=pval & stat_p_all(:,4)>=pval & stat_p_all(:,6)>=pval);
d_bsig = find(stat_p_all(:,3)<pval & stat_p_all(:,4)<pval);
d_esig = find(stat_p_all(:,3)<pval & stat_p_all(:,4)>=pval);
d_csig = find(stat_p_all(:,3)<pval & stat_p_all(:,4)>=pval & (stat_p_all(:,2)<pval | stat_p_all(:,6)<pval));
d_nsig = find(stat_p_all(:,3)<pval & stat_p_all(:,2)>=pval & stat_p_all(:,4)>=pval & stat_p_all(:,6)>=pval);
r_bsig = find(stat_p_all(:,5)<pval & stat_p_all(:,6)<pval);
r_esig = find(stat_p_all(:,5)<pval & stat_p_all(:,6)>=pval);
r_csig = find(stat_p_all(:,5)<pval & stat_p_all(:,6)>=pval & (stat_p_all(:,2)<pval | stat_p_all(:,4)<pval));
r_nsig = find(stat_p_all(:,5)<pval & stat_p_all(:,2)>=pval & stat_p_all(:,4)>=pval & stat_p_all(:,6)>=pval);

s_sig2 = find(stat_p_all(:,2)<pval);
d_sig2 = find(stat_p_all(:,4)<pval);
r_sig2 = find(stat_p_all(:,6)<pval);
s_bsig2 = find(stat_p_all(:,2)<pval & stat_p_all(:,1)<pval);
s_esig2 = find(stat_p_all(:,2)<pval & stat_p_all(:,1)>=pval);
s_csig2 = find(stat_p_all(:,2)<pval & stat_p_all(:,1)>=pval & (stat_p_all(:,3)<pval | stat_p_all(:,5)<pval));
s_nsig2 = find(stat_p_all(:,2)<pval & stat_p_all(:,1)>=pval & stat_p_all(:,3)>=pval & stat_p_all(:,5)>=pval);
d_bsig2 = find(stat_p_all(:,4)<pval & stat_p_all(:,3)<pval);
d_esig2 = find(stat_p_all(:,4)<pval & stat_p_all(:,3)>=pval);
d_csig2 = find(stat_p_all(:,4)<pval & stat_p_all(:,3)>=pval & (stat_p_all(:,1)<pval | stat_p_all(:,5)<pval));
d_nsig2 = find(stat_p_all(:,4)<pval & stat_p_all(:,1)>=pval & stat_p_all(:,3)>=pval & stat_p_all(:,5)>=pval);
r_bsig2 = find(stat_p_all(:,6)<pval & stat_p_all(:,5)<pval);
r_esig2 = find(stat_p_all(:,6)<pval & stat_p_all(:,5)>=pval);
r_csig2 = find(stat_p_all(:,6)<pval & stat_p_all(:,5)>=pval & (stat_p_all(:,1)<pval | stat_p_all(:,3)<pval));
r_nsig2 = find(stat_p_all(:,6)<pval & stat_p_all(:,1)>=pval & stat_p_all(:,3)>=pval & stat_p_all(:,5)>=pval);


clear pop_map pop_remap pop_sel_10 pop_sel_10_re

pop_map{1,1} = stat_p_all(s_sig,:);
pop_map{1,2} = stat_p_all(d_sig,:);
pop_map{1,3} = stat_p_all(r_sig,:);

pop_sel_s{1,1} = resp_all_b(s_sig,:);
pop_sel_s{1,2} = resp_all_s(s_sig,:);
pop_sel_s{1,3} = resp_all_d(s_sig,:);
pop_sel_s{1,4} = resp_all_r(s_sig,:);

pop_sel_d{1,1} = resp_all_b(d_sig,:);
pop_sel_d{1,2} = resp_all_s(d_sig,:);
pop_sel_d{1,3} = resp_all_d(d_sig,:);
pop_sel_d{1,4} = resp_all_r(d_sig,:);

pop_sel_r{1,1} = resp_all_b(r_sig,:);
pop_sel_r{1,2} = resp_all_s(r_sig,:);
pop_sel_r{1,3} = resp_all_d(r_sig,:);
pop_sel_r{1,4} = resp_all_r(r_sig,:);

pop_map{2,1} = stat_p_all(s_sig2,:);
pop_map{2,2} = stat_p_all(d_sig2,:);
pop_map{2,3} = stat_p_all(r_sig2,:);

pop_sel_s{2,1} = resp_all_b(s_sig2,:);
pop_sel_s{2,2} = resp_all_s(s_sig2,:);
pop_sel_s{2,3} = resp_all_d(s_sig2,:);
pop_sel_s{2,4} = resp_all_r(s_sig2,:);

pop_sel_d{2,1} = resp_all_b(d_sig2,:);
pop_sel_d{2,2} = resp_all_s(d_sig2,:);
pop_sel_d{2,3} = resp_all_d(d_sig2,:);
pop_sel_d{2,4} = resp_all_r(d_sig2,:);

pop_sel_r{2,1} = resp_all_b(r_sig2,:);
pop_sel_r{2,2} = resp_all_s(r_sig2,:);
pop_sel_r{2,3} = resp_all_d(r_sig2,:);
pop_sel_r{2,4} = resp_all_r(r_sig2,:);


clear temp aa index
temp = abs(resp_all_s(:,1));
[aa index] = sortrows(-temp);
pop_sel_10{1,1} = resp_all_s(index(1:round(length(index))/topcri),:);
clear temp aa index
temp = pop_sel_10{1,1}(:,1);
[aa index] = sortrows(-temp);
pop_sel_10_re{1,1} = pop_sel_10{1,1}(index,:);

clear temp aa index
temp = abs(resp_all_s(:,2));
[aa index] = sortrows(-temp);
pop_sel_10{1,2} = resp_all_s(index(1:round(length(index))/topcri),:);
clear temp aa index
temp = pop_sel_10{1,2}(:,2);
[aa index] = sortrows(-temp);
pop_sel_10_re{1,2} = pop_sel_10{1,2}(index,:);

clear temp aa index
temp = abs(resp_all_d(:,1));
[aa index] = sortrows(-temp);
pop_sel_10{2,1} = resp_all_d(index(1:round(length(index))/topcri),:);
clear temp aa index
temp = pop_sel_10{2,1}(:,1);
[aa index] = sortrows(-temp);
pop_sel_10_re{2,1} = pop_sel_10{2,1}(index,:);

clear temp aa index
temp = abs(resp_all_d(:,2));
[aa index] = sortrows(-temp);
pop_sel_10{2,2} = resp_all_d(index(1:round(length(index))/topcri),:);
clear temp aa index
temp = pop_sel_10{2,2}(:,2);
[aa index] = sortrows(-temp);
pop_sel_10_re{2,2} = pop_sel_10{2,2}(index,:);

clear temp aa index
temp = abs(resp_all_r(:,1));
[aa index] = sortrows(-temp);
pop_sel_10{3,1} = resp_all_r(index(1:round(length(index))/topcri),:);
clear temp aa index
temp = pop_sel_10{3,1}(:,1);
[aa index] = sortrows(-temp);
pop_sel_10_re{3,1} = pop_sel_10{3,1}(index,:);

clear temp aa index
temp = abs(resp_all_r(:,2));
[aa index] = sortrows(-temp);
pop_sel_10{3,2} = resp_all_r(index(1:round(length(index))/topcri),:);
clear temp aa index
temp = pop_sel_10{3,2}(:,2);
[aa index] = sortrows(-temp);
pop_sel_10_re{3,2} = pop_sel_10{3,2}(index,:);



for j=1:3
    clear temp1 temp2 aa bb index1 index2
    temp1 = pop_map{1,j}(:,2*j-1);
    temp2 = pop_map{2,j}(:,2*j);
    
    [aa index1] = sortrows(temp1);
    [bb index2] = sortrows(temp2);
    
    pop_remap{1,j} = pop_map{1,j}(index1,:);
    pop_remap{2,j} = pop_map{2,j}(index2,:);
    
end

clear pop_sel_s_re pop_sel_d_re pop_sel_r_re
for i=1:2
    for j=1:4
        clear temp1 aa1 index1
        clear temp2 aa2 index2
        clear temp3 aa3 index3
        temp1 = pop_sel_s{i,2}(:,i);
        temp2 = pop_sel_d{i,3}(:,i);
        temp3 = pop_sel_r{i,4}(:,i);
        
        [aa1 index1] = sortrows(-temp1);
        [aa2 index2] = sortrows(-temp2);
        [aa3 index3] = sortrows(-temp3);
        
        pop_sel_s_re{i,j} = pop_sel_s{i,j}(index1,:);
        pop_sel_d_re{i,j} = pop_sel_d{i,j}(index2,:);
        pop_sel_r_re{i,j} = pop_sel_r{i,j}(index3,:);
        
    end
end

% curr = 0;
% figure
% for i=1:2
%     for j=1:4
%         curr = curr + 1;
%         subplot(6,4,curr)
%         imagesc(pop_sel_s_re{i,j})
%         caxis([-1 1])
%         
%         subplot(6,4,curr+8)
%         scatter(pop_sel_s_re{i,j}(:,1),pop_sel_s_re{i,j}(:,2),sz,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.2)
%         hold on
%         title(strcat(num2str(length(pop_sel_s_re{i,j})),'/',num2str(length(resp_all_b))))
%         xlim([-1 1]) 
%         ylim([-1 1])
%         
%         subplot(6,4,curr+16)
%         clear aa bb cc
%         if i == 1
%             aa = histcounts(pop_sel_s_re{i,j}(:,1),edges1);
%             bb = histcounts(pop_sel_s_re{i,j}(find(pop_sel_s_re{i,j}(:,1)>0),2),edges1);
%             cc = histcounts(pop_sel_s_re{i,j}(find(pop_sel_s_re{i,j}(:,1)<0),2),edges1);
%         elseif i == 2
%             aa = histcounts(pop_sel_s_re{i,j}(:,2),edges1);
%             bb = histcounts(pop_sel_s_re{i,j}(find(pop_sel_s_re{i,j}(:,2)>0),1),edges1);
%             cc = histcounts(pop_sel_s_re{i,j}(find(pop_sel_s_re{i,j}(:,2)<0),1),edges1);
%         end
%         
%         hold on
%         bar(edges2,bb,'FaceColor','b','FaceAlpha',.5)
%         hold on
%         bar(edges2,cc,'FaceColor','r','FaceAlpha',.5)
%         hold on
%         plot(edges2,aa,'k','LineWidth',1)
%     end
% end
% 
% 
% curr = 0;
% figure
% for i=1:2
%     for j=1:4
%         curr = curr + 1;
%         subplot(6,4,curr)
%         imagesc(pop_sel_d_re{i,j})
%         caxis([-1 1])
%         
%         subplot(6,4,curr+8)
%         scatter(pop_sel_d_re{i,j}(:,1),pop_sel_d_re{i,j}(:,2),sz,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.2)
%         hold on
%         title(strcat(num2str(length(pop_sel_d_re{i,j})),'/',num2str(length(resp_all_b))))
%         xlim([-1 1]) 
%         ylim([-1 1])
%         
%         subplot(6,4,curr+16)
%         clear aa bb cc
%         if i == 1
%             aa = histcounts(pop_sel_d_re{i,j}(:,1),edges1);
%             bb = histcounts(pop_sel_d_re{i,j}(find(pop_sel_d_re{i,j}(:,1)>0),2),edges1);
%             cc = histcounts(pop_sel_d_re{i,j}(find(pop_sel_d_re{i,j}(:,1)<0),2),edges1);
%         elseif i == 2
%             aa = histcounts(pop_sel_d_re{i,j}(:,2),edges1);
%             bb = histcounts(pop_sel_d_re{i,j}(find(pop_sel_d_re{i,j}(:,2)>0),1),edges1);
%             cc = histcounts(pop_sel_d_re{i,j}(find(pop_sel_d_re{i,j}(:,2)<0),1),edges1);
%         end
%         
%         hold on
%         bar(edges2,bb,'FaceColor','b','FaceAlpha',.5)
%         hold on
%         bar(edges2,cc,'FaceColor','r','FaceAlpha',.5)
%         hold on
%         plot(edges2,aa,'k','LineWidth',1)
%     end
% end
% 
% curr = 0;
% figure
% for i=1:2
%     for j=1:4
%         curr = curr + 1;
%         subplot(6,4,curr)
%         imagesc(pop_sel_r_re{i,j})
%         caxis([-1 1])
%         
%         subplot(6,4,curr+8)
%         scatter(pop_sel_r_re{i,j}(:,1),pop_sel_r_re{i,j}(:,2),sz,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.2)
%         hold on
%         title(strcat(num2str(length(pop_sel_r_re{i,j})),'/',num2str(length(resp_all_b))))
%         xlim([-1 1]) 
%         ylim([-1 1])
%         
%         subplot(6,4,curr+16)
%         clear aa bb cc
%         if i == 1
%             aa = histcounts(pop_sel_r_re{i,j}(:,1),edges1);
%             bb = histcounts(pop_sel_r_re{i,j}(find(pop_sel_r_re{i,j}(:,1)>0),2),edges1);
%             cc = histcounts(pop_sel_r_re{i,j}(find(pop_sel_r_re{i,j}(:,1)<0),2),edges1);
%         elseif i == 2
%             aa = histcounts(pop_sel_r_re{i,j}(:,2),edges1);
%             bb = histcounts(pop_sel_r_re{i,j}(find(pop_sel_r_re{i,j}(:,2)>0),1),edges1);
%             cc = histcounts(pop_sel_r_re{i,j}(find(pop_sel_r_re{i,j}(:,2)<0),1),edges1);
%         end
%         
%         hold on
%         bar(edges2,bb,'FaceColor','b','FaceAlpha',.5)
%         hold on
%         bar(edges2,cc,'FaceColor','r','FaceAlpha',.5)
%         hold on
%         plot(edges2,aa,'k','LineWidth',1)
%     end
% end


figure
curr = 0;
for i=1:2
    for j=1:4        
        curr = curr + 1;
        if j == 2
        subplot(6,4,curr)
        imagesc(pop_sel_s_re{i,j})
        caxis([-1 1])
        
        subplot(6,4,curr+8)
        scatter(pop_sel_s_re{i,j}(:,1),pop_sel_s_re{i,j}(:,2),sz,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.1)
        hold on
        title(strcat(num2str(length(pop_sel_s_re{i,j})),'/',num2str(length(resp_all_b))))
        xlim([-1 1]) 
        ylim([-1 1])
        
        subplot(6,4,curr+16)
        clear aa bb cc
        if i == 1
            aa = histcounts(pop_sel_s_re{i,j}(:,1),edges1);
            bb = histcounts(pop_sel_s_re{i,j}(find(pop_sel_s_re{i,j}(:,1)>0),2),edges1);
            cc = histcounts(pop_sel_s_re{i,j}(find(pop_sel_s_re{i,j}(:,1)<0),2),edges1);
        elseif i == 2
            aa = histcounts(pop_sel_s_re{i,j}(:,2),edges1);
            bb = histcounts(pop_sel_s_re{i,j}(find(pop_sel_s_re{i,j}(:,2)>0),1),edges1);
            cc = histcounts(pop_sel_s_re{i,j}(find(pop_sel_s_re{i,j}(:,2)<0),1),edges1);
        end
        
        hold on
        bar(edges2,bb,'FaceColor','b','FaceAlpha',.5)
        hold on
        bar(edges2,cc,'FaceColor','r','FaceAlpha',.5)
        hold on
        plot(edges2,aa,'k','LineWidth',1)
        end
    end
end


curr = 0;
for i=1:2
    for j=1:4
        curr = curr + 1;
        if j == 3
        subplot(6,4,curr)
        imagesc(pop_sel_d_re{i,j})
        caxis([-1 1])
        
        subplot(6,4,curr+8)
        scatter(pop_sel_d_re{i,j}(:,1),pop_sel_d_re{i,j}(:,2),sz,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.1)
        hold on
        title(strcat(num2str(length(pop_sel_d_re{i,j})),'/',num2str(length(resp_all_b))))
        xlim([-1 1]) 
        ylim([-1 1])
        
        subplot(6,4,curr+16)
        clear aa bb cc
        if i == 1
            aa = histcounts(pop_sel_d_re{i,j}(:,1),edges1);
            bb = histcounts(pop_sel_d_re{i,j}(find(pop_sel_d_re{i,j}(:,1)>0),2),edges1);
            cc = histcounts(pop_sel_d_re{i,j}(find(pop_sel_d_re{i,j}(:,1)<0),2),edges1);
        elseif i == 2
            aa = histcounts(pop_sel_d_re{i,j}(:,2),edges1);
            bb = histcounts(pop_sel_d_re{i,j}(find(pop_sel_d_re{i,j}(:,2)>0),1),edges1);
            cc = histcounts(pop_sel_d_re{i,j}(find(pop_sel_d_re{i,j}(:,2)<0),1),edges1);
        end
        
        hold on
        bar(edges2,bb,'FaceColor','b','FaceAlpha',.5)
        hold on
        bar(edges2,cc,'FaceColor','r','FaceAlpha',.5)
        hold on
        plot(edges2,aa,'k','LineWidth',1)
        end
    end
end

curr = 0;
for i=1:2
    for j=1:4        
        curr = curr + 1;
        if j == 4
        subplot(6,4,curr)
        imagesc(pop_sel_r_re{i,j})
        caxis([-1 1])
        
        subplot(6,4,curr+8)
        scatter(pop_sel_r_re{i,j}(:,1),pop_sel_r_re{i,j}(:,2),sz,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.1)
        hold on
        title(strcat(num2str(length(pop_sel_r_re{i,j})),'/',num2str(length(resp_all_b))))
        xlim([-1 1]) 
        ylim([-1 1])
        
        subplot(6,4,curr+16)
        clear aa bb cc
        if i == 1
            aa = histcounts(pop_sel_r_re{i,j}(:,1),edges1);
            bb = histcounts(pop_sel_r_re{i,j}(find(pop_sel_r_re{i,j}(:,1)>0),2),edges1);
            cc = histcounts(pop_sel_r_re{i,j}(find(pop_sel_r_re{i,j}(:,1)<0),2),edges1);
        elseif i == 2
            aa = histcounts(pop_sel_r_re{i,j}(:,2),edges1);
            bb = histcounts(pop_sel_r_re{i,j}(find(pop_sel_r_re{i,j}(:,2)>0),1),edges1);
            cc = histcounts(pop_sel_r_re{i,j}(find(pop_sel_r_re{i,j}(:,2)<0),1),edges1);
        end
        
        hold on
        bar(edges2,bb,'FaceColor','b','FaceAlpha',.5)
        hold on
        bar(edges2,cc,'FaceColor','r','FaceAlpha',.5)
        hold on
        plot(edges2,aa,'k','LineWidth',1)
        end
    end
end


% figure
% curr = 0;
% for j=1:2
%     for i=1:3        
%         curr = curr + 1;
%         subplot(6,3,curr)
%         imagesc(pop_sel_10_re{i,j})
%         caxis([-1 1])
%         
%         subplot(6,3,curr+6)
%         scatter(pop_sel_10_re{i,j}(:,1),pop_sel_10_re{i,j}(:,2),sz,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.1)
%         hold on
%         title(strcat(num2str(length(pop_sel_10_re{i,j})),'/',num2str(length(stat_p_all))))
%         xlim([-1 1]) 
%         ylim([-1 1])
%         
%         subplot(6,3,curr+12)
%         clear aa bb cc
%         if j == 1
%             aa = histcounts(pop_sel_10_re{i,j}(:,1),edges1);
%             bb = histcounts(pop_sel_10_re{i,j}(find(pop_sel_10_re{i,j}(:,1)>0),2),edges1);
%             cc = histcounts(pop_sel_10_re{i,j}(find(pop_sel_10_re{i,j}(:,1)<0),2),edges1);
%         elseif j == 2
%             aa = histcounts(pop_sel_10_re{i,j}(:,2),edges1);
%             bb = histcounts(pop_sel_10_re{i,j}(pop_sel_10_re{i,j}(:,2)>0,1),edges1);
%             cc = histcounts(pop_sel_10_re{i,j}(find(pop_sel_10_re{i,j}(:,2)<0),1),edges1);
%         end
%         
%         hold on
%         bar(edges2,bb,'FaceColor','b','FaceAlpha',.5)
%         hold on
%         bar(edges2,cc,'FaceColor','r','FaceAlpha',.5)
%         hold on
%         plot(edges2,aa,'k','LineWidth',1)
%     end
% end

clear pop_sel_stab
pop_sel_stab{1,1} = pop_sel_s_re{1,2};
pop_sel_stab{2,1} = pop_sel_s_re{2,2};
pop_sel_stab{1,2} = pop_sel_d_re{1,3};
pop_sel_stab{2,2} = pop_sel_d_re{2,3};
pop_sel_stab{1,3} = pop_sel_r_re{1,4};
pop_sel_stab{2,3} = pop_sel_r_re{2,4};

ti_text = ["Sample","Delay","Response"];

edges1 = [-1.05:0.1:1.05];
edges2 = [-1:0.1:1];

% segmentations depending on the amplitude of selectivity index
clear mean_si std_si
clear seg_text
seg_text = {};

for t=1:segno
    seg_text{t} = [num2str(t),'/',num2str(segno)];
    
end

figure
for j=1:3
    clear temp1 aa index segcri temp2
    temp1 = pop_sel_stab{1,j};
    temp1(isnan(temp1))=0;
    [aa index] = sortrows(-abs(temp1(:,1)));
    
    segcri = [1:floor(length(aa)/segno)-1:length(aa)];
    
    for t=1:segno
        temp2 = temp1(index(segcri(t):segcri(t+1)-1),:);
        clear aaa bbb ccc
        aaa = histcounts(temp2(:,1),edges1);
        bbb = histcounts(temp2(find(temp2(:,1)>0),2),edges1);
        ccc = histcounts(temp2(find(temp2(:,1)<0),2),edges1);
        
        subplot(3,segno+1,(segno+1)*(j-1)+t)
        hold on
        bar(edges2,bbb,'FaceColor','b','FaceAlpha',.6)
        hold on
        bar(edges2,ccc,'FaceColor','r','FaceAlpha',.6)
        hold on
        plot(edges2,aaa,'k','LineWidth',1.5)
        hold on
        title(strcat(ti_text(j),'-',num2str(t),'/',num2str(segno)))
        xlabel('selectivity index')
        ylabel('cell no.')
        
        mean_si{j,1}(1,t) = mean(temp2(find(temp2(:,1)>0),2));
        std_si{j,1}(1,t) = std(temp2(find(temp2(:,1)>0),2))/sqrt(length(bbb));
        mean_si{j,1}(2,t) = mean(temp2(find(temp2(:,1)<0),2));
        std_si{j,1}(2,t) = std(temp2(find(temp2(:,1)>0),2))/sqrt(length(ccc));
        
        mean_si{j,1}(3,t) = mean(temp2(find(temp2(:,1)>0),1));
        std_si{j,1}(3,t) = std(temp2(find(temp2(:,1)>0),1))/sqrt(length(bbb));
        mean_si{j,1}(4,t) = mean(temp2(find(temp2(:,1)<0),1));
        std_si{j,1}(4,t) = std(temp2(find(temp2(:,1)>0),1))/sqrt(length(ccc));
    end
    
    subplot(3,segno+1,(segno+1)*(j-1)+segno+1)
    hold on
    errorbar(mean_si{j,1}(1,:),std_si{j,1}(1,:),'b')
    hold on
    errorbar(mean_si{j,1}(2,:),std_si{j,1}(2,:),'r')
    hold on
    errorbar(mean_si{j,1}(3,:),std_si{j,1}(3,:),'k')
    hold on
    errorbar(mean_si{j,1}(4,:),std_si{j,1}(4,:),'k')
    ylim([-1 1])
    xlim([0 segno+1])
    xticks([1:1:segno]);
    xticklabels(seg_text)
    xlabel('segments')
    ylabel('mean SI')
end
hold on
sgtitle('Early -> Late')

figure
for j=1:3
    clear temp1 aa index segcri temp2
    temp1 = pop_sel_stab{2,j};
    temp1(isnan(temp1))=0;
    [aa index] = sortrows(-abs(temp1(:,2)));
    
    segcri = [1:floor(length(aa)/segno)-1:length(aa)];
    
    for t=1:segno
        temp2 = temp1(index(segcri(t):segcri(t+1)-1),:);
        clear aaa bbb ccc
        aaa = histcounts(temp2(:,2),edges1);
        bbb = histcounts(temp2(find(temp2(:,2)>0),1),edges1);
        ccc = histcounts(temp2(find(temp2(:,2)<0),1),edges1);
        
        subplot(3,segno+1,(segno+1)*(j-1)+t)
        hold on
        bar(edges2,bbb,'FaceColor','b','FaceAlpha',.6)
        hold on
        bar(edges2,ccc,'FaceColor','r','FaceAlpha',.6)
        hold on
        plot(edges2,aaa,'k','LineWidth',1.5)
        hold on
        title(strcat(ti_text(j),'-',num2str(t),'/',num2str(segno)))
        xlabel('selectivity index')
        ylabel('cell no.')
        
        mean_si{j,2}(1,t) = mean(temp2(find(temp2(:,2)>0),1));
        std_si{j,2}(1,t) = std(temp2(find(temp2(:,2)>0),1))/sqrt(length(bbb));
        mean_si{j,2}(2,t) = mean(temp2(find(temp2(:,2)<0),1));
        std_si{j,2}(2,t) = std(temp2(find(temp2(:,2)>0),1))/sqrt(length(ccc));        
        
        mean_si{j,2}(3,t) = mean(temp2(find(temp2(:,2)>0),2));
        std_si{j,2}(3,t) = std(temp2(find(temp2(:,2)>0),2))/sqrt(length(bbb));
        mean_si{j,2}(4,t) = mean(temp2(find(temp2(:,2)<0),2));
        std_si{j,2}(4,t) = std(temp2(find(temp2(:,2)>0),2))/sqrt(length(ccc));
    end
    
    subplot(3,segno+1,(segno+1)*(j-1)+segno+1)
    hold on
    errorbar(mean_si{j,2}(1,:),std_si{j,2}(1,:),'b')
    hold on
    errorbar(mean_si{j,2}(2,:),std_si{j,2}(2,:),'r')
    hold on
    errorbar(mean_si{j,2}(3,:),std_si{j,2}(3,:),'k')
    hold on
    errorbar(mean_si{j,2}(4,:),std_si{j,2}(4,:),'k')
    ylim([-1 1])
    xlim([0 segno+1])
    xticks([1:1:segno]);
    xticklabels(seg_text)
    xlabel('segments')
    ylabel('mean SI')
end       
hold on
sgtitle('Late -> Early')  


figure
for j=1:3
    clear temp1 aa index segcri temp2    
    temp1 = pop_sel_stab{1,j};    
    temp1(isnan(temp1))=0;
    [aa index] = sortrows(-abs(temp1(:,1)));
    
    clear aa dd bb cc
    aa = histcounts(temp1(find(temp1(:,1)>0),1),edges1);
    dd = histcounts(temp1(find(temp1(:,1)<0),1),edges1);
    bb = histcounts(temp1(find(temp1(:,1)>0),2),edges1);
    cc = histcounts(temp1(find(temp1(:,1)<0),2),edges1);

    subplot(3,2,(j-1)*2+1)
    hold on
    bar(edges2,aa,'FaceColor','b','FaceAlpha',.9)
    hold on
    bar(edges2,dd,'FaceColor','r','FaceAlpha',.9)
    xlabel('selectivity index')
    ylabel('cell no.')
    ax=gca;
    ax.XAxis.FontSize=13;
    ax.YAxis.FontSize=13;
        
    subplot(3,2,j*2)
    hold on
    bar(edges2,bb,'FaceColor','b','FaceAlpha',.3)
    hold on
    bar(edges2,cc,'FaceColor','r','FaceAlpha',.3)
    xlabel('selectivity index')
    ylabel('cell no.')
    ax=gca;
    ax.XAxis.FontSize=13;
    ax.YAxis.FontSize=13;
    
end

hold on
sgtitle('Early -> Late')  

set(gcf,'color','w')

clear FN cri
cri = num2str(pval);
FN = strcat('within_context_SI_EL_p',cri(3:end),'.emf');
exportgraphics(gcf,FN,'ContentType','vector')

figure
for j=1:3
    clear temp1 aa index segcri temp2    
    temp1 = pop_sel_stab{2,j};    
    temp1(isnan(temp1))=0;
    [aa index] = sortrows(-abs(temp1(:,2)));
    
    clear aa dd bb cc
    aa = histcounts(temp1(find(temp1(:,2)>0),1),edges1);
    dd = histcounts(temp1(find(temp1(:,2)<0),1),edges1);
    bb = histcounts(temp1(find(temp1(:,2)>0),2),edges1);
    cc = histcounts(temp1(find(temp1(:,2)<0),2),edges1);

    subplot(3,2,(j-1)*2+1)
    hold on
    bar(edges2,aa,'FaceColor','b','FaceAlpha',.3)
    hold on
    bar(edges2,dd,'FaceColor','r','FaceAlpha',.3)
    xlabel('selectivity index')
    ylabel('cell no.')
    ax=gca;
    ax.XAxis.FontSize=13;
    ax.YAxis.FontSize=13;
        
    subplot(3,2,j*2)
    hold on
    bar(edges2,bb,'FaceColor','b','FaceAlpha',.9)
    hold on
    bar(edges2,cc,'FaceColor','r','FaceAlpha',.9)
    xlabel('selectivity index')
    ylabel('cell no.')
    ax=gca;
    ax.XAxis.FontSize=13;
    ax.YAxis.FontSize=13;
    
end

hold on
sgtitle('Late -> Early')  
   
set(gcf,'color','w')

clear FN cri
cri = num2str(pval);
FN = strcat('within_context_SI_LE_p',cri(3:end),'.emf');
exportgraphics(gcf,FN,'ContentType','vector')
%%
%exportgraphics(gcf,'across_contexts_SI_EL_p01.emf','ContentType','vector')

        
    %%

comp{1,1} = find(stat_p_all(s_sig,2)<pval);
comp{1,2} = find(stat_p_all(s_sig,4)<pval);
comp{1,3} = find(stat_p_all(s_sig,6)<pval);

comp{1,4} = s_sig;
comp{1,5} = s_bsig;
comp{1,6} = s_esig;
comp{1,7} = s_csig;
comp{1,8} = s_nsig;

comp{2,1} = find(stat_p_all(d_sig,2)<pval);
comp{2,2} = find(stat_p_all(d_sig,4)<pval);
comp{2,3} = find(stat_p_all(d_sig,6)<pval);

comp{2,4} = d_sig;
comp{2,5} = d_bsig;
comp{2,6} = d_esig;
comp{2,7} = d_csig;
comp{2,8} = d_nsig;

comp{3,1} = find(stat_p_all(r_sig,2)<pval);
comp{3,2} = find(stat_p_all(r_sig,4)<pval);
comp{3,3} = find(stat_p_all(r_sig,6)<pval);

comp{3,4} = r_sig;
comp{3,5} = r_bsig;
comp{3,6} = r_esig;
comp{3,7} = r_csig;
comp{3,8} = r_nsig;

comp2{1,1} = find(stat_p_all(s_sig2,1)<pval);
comp2{1,2} = find(stat_p_all(s_sig2,3)<pval);
comp2{1,3} = find(stat_p_all(s_sig2,5)<pval);

comp2{1,4} = s_sig2;
comp2{1,5} = s_bsig2;
comp2{1,6} = s_esig2;
comp2{1,7} = s_csig2;
comp2{1,8} = s_nsig2;

comp2{2,1} = find(stat_p_all(d_sig2,1)<pval);
comp2{2,2} = find(stat_p_all(d_sig2,3)<pval);
comp2{2,3} = find(stat_p_all(d_sig2,5)<pval);

comp2{2,4} = d_sig2;
comp2{2,5} = d_bsig2;
comp2{2,6} = d_esig2;
comp2{2,7} = d_csig2;
comp2{2,8} = d_nsig2;

comp2{3,1} = find(stat_p_all(r_sig2,1)<pval);
comp2{3,2} = find(stat_p_all(r_sig2,3)<pval);
comp2{3,3} = find(stat_p_all(r_sig2,5)<pval);

comp2{3,4} = r_sig2;
comp2{3,5} = r_bsig2;
comp2{3,6} = r_esig2;
comp2{3,7} = r_csig2;
comp2{3,8} = r_nsig2;

for i=1:size(comp2,1)
    for j=1:size(comp2,2)
        takeaway1(i,j)=length(comp{i,j});
        takeaway2(i,j)=length(comp2{i,j});
    end
end




