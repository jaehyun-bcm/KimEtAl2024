
save responsiveness.mat resp_all_b resp_all_s resp_all_d resp_all_r


%%
clear
load responsiveness.mat
%%
sz = 1.5; 

clear colorcri
% across contexts
maxcri = 80;
mincri = 13;
% within context over time
%maxcri = 50;
%mincri = 8;
colorcri = parula(maxcri-mincri);
edges1 = [-1:.05:1];

clear targetdata metrics_matrix_s metrics_matrix_d metrics_matrix_r
targetdata1 = resp_all_s;
targetdata2 = resp_all_d;
targetdata3 = resp_all_r;

for i=1:size(edges1,2)-1         
    for j=1:size(edges1,2)-1
        clear currpx
        currpx = length(find(targetdata1(:,1) > edges1(i) & targetdata1(:,1) <= edges1(i+1) & targetdata1(:,2) > edges1(j) & targetdata1(:,2) <= edges1(j+1)))- mincri;
        if currpx > maxcri-mincri
            currpx = maxcri-mincri;
        end
        
        if currpx < 1
            metrics_matrix_s(j,i,1:3) = [1 1 1];
        else 
            metrics_matrix_s(j,i,1:3) = colorcri(currpx,:);
        end
        
        clear currpx
        currpx = length(find(targetdata2(:,1) > edges1(i) & targetdata2(:,1) <= edges1(i+1) & targetdata2(:,2) > edges1(j) & targetdata2(:,2) <= edges1(j+1)))- mincri;
        if currpx > maxcri-mincri
            currpx = maxcri-mincri;
        end
        
        if currpx < 1
            metrics_matrix_d(j,i,1:3) = [1 1 1];
        else 
            metrics_matrix_d(j,i,1:3) = colorcri(currpx,:);
        end
        
        clear currpx
        currpx = length(find(targetdata3(:,1) > edges1(i) & targetdata3(:,1) <= edges1(i+1) & targetdata3(:,2) > edges1(j) & targetdata3(:,2) <= edges1(j+1)))- mincri;
        if currpx > maxcri-mincri
            currpx = maxcri-mincri;
        end
        
        if currpx < 1
            metrics_matrix_r(j,i,1:3) = [1 1 1];
        else 
            metrics_matrix_r(j,i,1:3) = colorcri(currpx,:);
        end
    end
end

figure
subplot(2,3,1)
hold on
scatter(resp_all_s(:,1),resp_all_s(:,2),sz,'k')

subplot(2,3,2)
hold on
scatter(resp_all_d(:,1),resp_all_d(:,2),sz,'k')

subplot(2,3,3)
hold on
scatter(resp_all_r(:,1),resp_all_r(:,2),sz,'k')


subplot(2,3,4)
hold on
imagesc(metrics_matrix_s)

subplot(2,3,5)
hold on
imagesc(metrics_matrix_d)

subplot(2,3,6)
hold on
imagesc(metrics_matrix_r)

for i=1:3
    subplot(2,3,i)
    hold on
    line([-1 1],[0 0],'color','k')
    line([0 0],[-1 1],'color','k')
    xlim([-1 1])
    ylim([-1 1])
    xticks([-1:1:1])
    yticks([-1:1:1])
    xlabel('Weight-Context 1','fontsize',12)
    ylabel('Weight-Context 2','fontsize',12)
    %xlabel('W_E_a_r_l_y','fontsize',12)
    %ylabel('W_L_a_t_e','fontsize',12)
    
    subplot(2,3,3+i)
    line([0 length(edges1)],[length(edges1)/2 length(edges1)/2],'color','k')
    line([length(edges1)/2 length(edges1)/2],[0 length(edges1)],'color','k')
    xticks([0:(length(edges1)-1)/2:length(edges1)])
    xticklabels([-1:1:1])
    yticks([0:(length(edges1)-1)/2:length(edges1)])
    yticklabels([-1:1:1])
    xlim([0 length(edges1)])
    ylim([0 length(edges1)])
end

ax=gca;
ax.XAxis.FontSize=11;
ax.YAxis.FontSize=11;

sgtitle(strcat('maxcri:',num2str(maxcri),'mincri:',num2str(mincri),' - bins:',num2str(length(edges1)),' - cells: n=',num2str(size(resp_all_r,1))))
%%
set(gcf,'color','w')
exportgraphics(gcf,'acrosscontexts_SelectivityIndexScatter.emf','ContentType','vector')