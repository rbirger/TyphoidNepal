%% Data on the serovar (typhi vs paratyphi) is missing for June 2001-May 2002, so let's impute it based on the proportion typhi vs paratyphi using the rest of the data

typhi_nepal=C_nepal(:,3);
paratyphi_nepal=C_nepal(:,2);

nepal_ptyphi=sum(C_nepal([1:216 264:end],3))/sum(C_nepal([1:216 264:end],1));
nepal_pparatyphi=sum(C_nepal([1:216 264:end],2))/sum(C_nepal([1:216 264:end],1));

for i=217:263
    typhi_nepal(i,1)=poissrnd(nepal_ptyphi*C_nepal(i,1));
    paratyphi_nepal(i,1)=C_nepal(i,1)-typhi_nepal(i,1);
end

%%
figure
hold on
plot(datenum(date_nepal),typhi_nepal,'r')
%plot(datenum(date_nepal),C_nepal(:,3))
%plot(datenum(date_nepal),C_nepal(:,2),'Color',[0 .5 0])
datetick('x','mmm-yy')
xlim([datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7])
set(gca,'FontSize',14)
ylabel('Number of TF cases (per week)','FontSize',14)

%%
figure
hold on
plot(datenum(date_nepal(1:216,:)),C_nepal(1:216,3),'b')
plot(datenum(date_nepal(1:216,:)),C_nepal(1:216,2),'Color',[0 .5 0])
plot(datenum(date_nepal(217:263,:)),round(C_nepal(217:263,1)*nepal_ptyphi),'b') %,'--b')
plot(datenum(date_nepal(217:263,:)),round(C_nepal(217:263,1)*nepal_pparatyphi),'Color',[0 .5 0]) %,'--','Color',[0 .5 0])
plot(datenum(date_nepal(264:end,:)),C_nepal(264:end,3),'b')
plot(datenum(date_nepal(264:end,:)),C_nepal(264:end,2),'Color',[0 .5 0])
datetick('x','yyyy')
xlim([datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7])
ylabel('No. of cases (per week)')
legend('\it{S}.\rm Typhi','\it{S}.\rm Paratyphi')

%%
figure
subplot(3,1,1)
plot(datenum(date_nepal),typhi_nepal,'b')
datetick('x','mmm-yy')
xlim([datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7])
set(gca,'YLim',[0 80],'YTick',0:20:80)
ylabel({'No. of Typhi cases';'(per week)'})

subplot(3,1,2)
plot(datenum(date_nepal),paratyphi_nepal,'Color',[0 .5 0])
datetick('x','mmm-yy')
xlim([datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7])
set(gca,'YLim',[0 80],'YTick',0:20:80)
ylabel({'No. of Paratyphi cases';'(per week)'})

subplot(3,1,3)
plot(datenum(date_rainwk),rainfall_nepalwk,'c')
datetick('x','mmm-yy')
xlim([datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7])
ylabel('Rainfall (mm/week)')

%%
figure
[ax,y1,y2]=plotyy(datenum(date_nepal),typhi_nepal,datenum(date_rainwk),rainfall_nepalwk);
datetick('x','mmm-yy')
set(y1,'Color','r')
set(y2,'Color','b')
set(ax(1),'YColor','k','XLim',[datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7],'YTick',0:20:100)
set(ax(2),'YColor','k','XLim',[datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7],'XTick',[],'YTick',0:100:500)
ylabel('Typhoid cases (per week)')
set(get(ax(2),'YLabel'),'String','Rainfall (mm/wk)')
legend('Typhoid cases','Rainfall')

%%
figure
subplot(2,2,1)
[ax,y1,y2]=plotyy(datenum(date_nepal),typhi_nepal,datenum(date_rainwk),rainfall_nepalwk);
datetick('x','mmm-yy')
set(y1,'Color','r')
set(y2,'Color','b')
set(ax(1),'YColor','k','XLim',[datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7],'YTick',0:20:100)
set(ax(2),'YColor','k','XLim',[datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7],'XTick',[],'YTick',0:100:500)
ylabel('Typhoid cases (per week)')
set(get(ax(2),'YLabel'),'String','Rainfall (mm/wk)')
legend('Typhoid cases','Rainfall')

subplot(2,2,2)
[ax,y1,y2]=plotyy(datenum(date_nepal),paratyphi_nepal,datenum(date_rainwk),rainfall_nepalwk);
datetick('x','mmm-yy')
set(y1,'Color','r')
set(y2,'Color','b')
set(ax(1),'YColor','k','XLim',[datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7],'YTick',0:20:100)
set(ax(2),'YColor','k','XLim',[datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7],'XTick',[],'YTick',0:100:500)
ylabel('Paratyphoid cases (per week)')
set(get(ax(2),'YLabel'),'String','Rainfall (mm/wk)')
legend('Paratyphoid cases','Rainfall')

subplot(2,2,3)
scatter(rainfall_nepalwk,typhi_nepal(1:length(rainfall_nepalwk),1))
xlabel('Rainfall (mm/wk)')
ylabel('Typhoid cases (per week)')

subplot(2,2,4)
scatter(rainfall_nepalwk,paratyphi_nepal(1:length(rainfall_nepalwk),1))
xlabel('Rainfall (mm/wk)')
ylabel('Paratyphoid cases (per week)')


%%
figure
%[ax,y1,y2]=plotyy(1:52,mean(typhi_nepal_bywk,2),1:52,rainbywk_nepal);
[ax,y1,y2]=plotyy(1:52,mean(typhi_wkyr,2),1:52,rainbywk_nepal);
set(y1,'Color','r')
set(y2,'Color','b')
set(ax(1),'YColor','k','XLim',[0 53],'XTick',0:13:52)
set(ax(2),'YColor','k','XLim',[0 53],'XTick',0:13:52)
ylabel('Avg typhoid cases (per week)')
set(get(ax(2),'YLabel'),'String','Avg rainfall (mm/wk)')
legend('Typhoid cases','Rainfall')


%%
figure
%[ax,y1,y2]=plotyy(1:52,mean(typhi_nepal_bywk,2),1:52,rainbywk_nepal);
[ax,y1,y2]=plotyy(1:52,typhi_wkyr,1:52,rainbywk_nepal);
set(y1,'Color','r')
set(y2,'Color','b')
set(ax(1),'YColor','k','XLim',[0 53],'XTick',0:13:52)
set(ax(2),'YColor','k','XLim',[0 53],'XTick',0:13:52)
ylabel('Avg typhoid cases (per week)')
set(get(ax(2),'YLabel'),'String','Avg rainfall (mm/wk)')
legend('Typhoid cases','Rainfall')