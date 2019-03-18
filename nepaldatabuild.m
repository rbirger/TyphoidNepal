%% Data on the serovar (typhi vs paratyphi) is missing for June 2001-May 2002, so let's impute it based on the proportion typhi vs paratyphi using the rest of the data
dbstop if error
typhi_nepal=C_nepal(:,3);
paratyphi_nepal=C_nepal(:,2);

nepal_ptyphi=sum(C_nepal([1:216 264:end],3))/sum(C_nepal([1:216 264:end],1));
nepal_pparatyphi=sum(C_nepal([1:216 264:end],2))/sum(C_nepal([1:216 264:end],1));

for i=217:263
    typhi_nepal(i,1)=poissrnd(nepal_ptyphi*C_nepal(i,1));
    paratyphi_nepal(i,1)=max(C_nepal(i,1)-typhi_nepal(i,1),0);
end

% First, log-transform the data (try also --and scale it to have a mean=0
% and variance=1)
logtyphi=log(typhi_nepal+1); 
lograinfall=log(rainfall_nepalwk+1);

%% Calculate correlation between rainfall and typhi/paratyphi cases and different lags (0-8 weeks)

corr_raincases=zeros(9,2);

for i=0:8
    corr_raincases(i+1,1)=corr(rainfall_nepalwk,typhi_nepal(1+i:length(rainfall_nepalwk)+i,1));
    corr_raincases(i+1,2)=corr(rainfall_nepalwk,paratyphi_nepal(1+i:length(rainfall_nepalwk)+i,1));
end

%% Fit harmonic regression

time=(15:length(date_nepal)+14)'; % data starts April 13 = week 15

mdl_typhi1=GeneralizedLinearModel.fit([time cos(2*pi*time/52.18) sin(2*pi*time/52.18)],logtyphi,'linear','Distribution','poisson','Link','log');
mdl_typhi2=GeneralizedLinearModel.fit([time cos(2*pi*time/52.18) sin(2*pi*time/52.18) cos(2*pi*time/26.09) sin(2*pi*time/26.09)],logtyphi,'linear','Distribution','poisson','Link','log');
%mdl_typhi3=GeneralizedLinearModel.fit([time cos(2*pi*time/52.18) sin(2*pi*time/52.18) cos(2*pi*time/26.09) sin(2*pi*time/26.09) cos(2*pi*time/13) sin(2*pi*time/13)],typhi_nepal,'linear','Distribution','poisson','Link','log');
%mdl_typhi4=GeneralizedLinearModel.fit([time cos(2*pi*time/156.5) sin(2*pi*time/156.5) cos(2*pi*time/52.18) sin(2*pi*time/52.18) cos(2*pi*time/26.09) sin(2*pi*time/26.09) cos(2*pi*time/13) sin(2*pi*time/13)],typhi_nepal,'linear','Distribution','poisson','Link','log');
%mdl_typhi5=GeneralizedLinearModel.fit([time cos(2*pi*time/156.5) sin(2*pi*time/156.5) cos(2*pi*time/52.18) sin(2*pi*time/52.18) cos(2*pi*time/26.09) sin(2*pi*time/26.09)],typhi_nepal,'linear','Distribution','poisson','Link','log');
%mdl_typhi6=GeneralizedLinearModel.fit([time cos(2*pi*time/521.8) sin(2*pi*time/521.8) cos(2*pi*time/156.5) sin(2*pi*time/156.5) cos(2*pi*time/52.18) sin(2*pi*time/52.18) cos(2*pi*time/26.09) sin(2*pi*time/26.09)],typhi_nepal,'linear','Distribution','poisson','Link','log');
mdl_typhi7=GeneralizedLinearModel.fit([time cos(2*pi*time/626.1) sin(2*pi*time/626.1) cos(2*pi*time/156.5) sin(2*pi*time/156.5) cos(2*pi*time/52.18) sin(2*pi*time/52.18) cos(2*pi*time/26.09) sin(2*pi*time/26.09)],logtyphi,'linear','Distribution','poisson','Link','log');
%mdl_typhi8=GeneralizedLinearModel.fit([time cos(2*pi*time/(52.18*14)) sin(2*pi*time/(52.18*14)) cos(2*pi*time/156.5) sin(2*pi*time/156.5) cos(2*pi*time/52.18) sin(2*pi*time/52.18) cos(2*pi*time/26.09) sin(2*pi*time/26.09)],typhi_nepal,'linear','Distribution','poisson','Link','log');

mdl_typhi3=GeneralizedLinearModel.fit([time cos(2*pi*time/156.5) sin(2*pi*time/156.5) cos(2*pi*time/52.18) sin(2*pi*time/52.18)],logtyphi,'linear','Distribution','poisson','Link','log');
mdl_typhi4=GeneralizedLinearModel.fit([time cos(2*pi*time/626.1) sin(2*pi*time/626.1) cos(2*pi*time/52.18) sin(2*pi*time/52.18)],logtyphi,'linear','Distribution','poisson','Link','log');
mdl_typhi5=GeneralizedLinearModel.fit([time cos(2*pi*time/626.1) sin(2*pi*time/626.1) cos(2*pi*time/156.5) sin(2*pi*time/156.5) cos(2*pi*time/52.18) sin(2*pi*time/52.18)],logtyphi,'linear','Distribution','poisson','Link','log');

%% display AICs
AICs_typhoid = [];
for i = [1:5,7]
   eval(['tmpMdl = mdl_typhi',num2str(i),';'])
   %disp(i);
   %disp(tmpMdl.ModelCriterion.AIC);
   AICs_typhoid = [AICs_typhoid;tmpMdl.ModelCriterion.AIC];
end


%%
%model 4 has lowest AIC
%Y_typhi6=predict(mdl_typhi6);
Y_typhi7=predict(mdl_typhi7);
%Y_typhi8=predict(mdl_typhi8);
Y_typhi4=max(0,exp(predict(mdl_typhi4))-1);

figure; hold on
plot(datenum(date_nepal),typhi_nepal,'b')
plot(datenum(date_nepal),Y_typhi4,'r')
%plot(datenum(date_nepal),Y_typhi6,'r')
plot(datenum(date_nepal),Y_typhi7,'g')
%plot(datenum(date_nepal),Y_typhi8,'c')
datetick('x','mmm-yy')
xlim([datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7])
set(gca,'YLim',[0 80],'YTick',0:20:80)
ylabel({'No. of Typhi cases';'(per week)'})
legend('Typhoid Cases','Model Fit 4 (12, 1 year periods)','Model Fit 7 (12,3,1, .5 yr periods)')

%%
% change residuals so they are divided, for log
%resid_typhi=typhi_nepal-Y_typhi4;

resid_typhi=typhi_nepal./Y_typhi4;

figure
[ax,y1,y2]=plotyy(datenum(date_nepal),resid_typhi,datenum(date_rainwk),rainfall_nepalwk);
datetick('x','mmm-yy')
set(y1,'Color','r')
set(y2,'Color','b')
set(ax(1),'YColor','k','XLim',[datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7],'YTick',0:2:10)
set(ax(2),'YColor','k','XLim',[datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7],'XTick',[],'YTick',0:100:500)
ylabel('Typhoid Residuals (per week)')
set(get(ax(2),'YLabel'),'String','Rainfall (mm/wk)')
legend('Typhoid Residuals','Rainfall')

%% Calculate correlation between rainfall and typhi/paratyphi cases and different lags (0-8 weeks)

corr_rainresid=zeros(9,1);

for i=0:8
    corr_rainresid(i+1,1)=corr(rainfall_nepalwk,resid_typhi(1+i:length(rainfall_nepalwk)+i,1));
    %corr_rainresid(i+1,2)=corr(rainfall_nepalwk,resid_paratyphi_nepal(1+i:length(rainfall_nepalwk)+i,1));
end

%% Fit harmonic regression to rainfall data 
timerain=time(1:716);

mdl_rainfall1=GeneralizedLinearModel.fit([timerain cos(2*pi*timerain/52.18) sin(2*pi*timerain/52.18)],rainfall_nepalwk,'linear','Distribution','poisson','Link','log');
mdl_rainfall2=GeneralizedLinearModel.fit([timerain cos(2*pi*timerain/52.18) sin(2*pi*timerain/52.18) cos(2*pi*timerain/26.09) sin(2*pi*timerain/26.09)],rainfall_nepalwk,'linear','Distribution','poisson','Link','log');
mdl_rainfall3=GeneralizedLinearModel.fit([timerain cos(2*pi*timerain/52.18) sin(2*pi*timerain/52.18) cos(2*pi*timerain/26.09) sin(2*pi*timerain/26.09) cos(2*pi*timerain/13) sin(2*pi*timerain/13)],rainfall_nepalwk,'linear','Distribution','poisson','Link','log');
mdl_rainfall4=GeneralizedLinearModel.fit([timerain cos(2*pi*timerain/156.5) sin(2*pi*timerain/156.5) cos(2*pi*timerain/52.18) sin(2*pi*timerain/52.18) cos(2*pi*timerain/26.09) sin(2*pi*timerain/26.09) cos(2*pi*timerain/13) sin(2*pi*timerain/13)],rainfall_nepalwk,'linear','Distribution','poisson','Link','log');
mdl_rainfall5=GeneralizedLinearModel.fit([timerain cos(2*pi*timerain/156.5) sin(2*pi*timerain/156.5) cos(2*pi*timerain/52.18) sin(2*pi*timerain/52.18) cos(2*pi*timerain/26.09) sin(2*pi*timerain/26.09)],rainfall_nepalwk,'linear','Distribution','poisson','Link','log');
mdl_rainfall6=GeneralizedLinearModel.fit([timerain cos(2*pi*timerain/521.8) sin(2*pi*timerain/521.8) cos(2*pi*timerain/156.5) sin(2*pi*timerain/156.5) cos(2*pi*timerain/52.18) sin(2*pi*timerain/52.18) cos(2*pi*timerain/26.09) sin(2*pi*timerain/26.09)],rainfall_nepalwk,'linear','Distribution','poisson','Link','log');
mdl_rainfall7=GeneralizedLinearModel.fit([timerain cos(2*pi*timerain/626.1) sin(2*pi*timerain/626.1) cos(2*pi*timerain/156.5) sin(2*pi*timerain/156.5) cos(2*pi*timerain/52.18) sin(2*pi*timerain/52.18) cos(2*pi*timerain/26.09) sin(2*pi*timerain/26.09)],rainfall_nepalwk,'linear','Distribution','poisson','Link','log');

%% Compare AICs for rainfall
AICs_rain = [];
for i = [1:7]
   eval(['tmpMdl = mdl_rainfall',num2str(i),';'])
   %disp(i);
   %disp(tmpMdl.ModelCriterion.AIC);
   AICs_rain = [AICs_rain;tmpMdl.ModelCriterion.AIC];
end

%%
%change to just weekly
%Y_rainfall7=predict(mdl_rainfall7);
%resid_rainfall=rainfall_nepalwk-Y_rainfall7;

% change to best model, #6
% and change to dividing for residuals instead of subtracting

Y_rainfall6=predict(mdl_rainfall6);
resid_rainfall=rainfall_nepalwk./Y_rainfall6;

% define lag index vector
lag_ind= [1:716];
date_num_nepal = datenum(date_nepal);
figure; 
subplot(1,2,1); hold on;
plot(datenum(date_rainwk),resid_rainfall,'b')
%plot(datenum(date_nepal),resid_typhi,'r')
plot(date_num_nepal(lag_ind),resid_typhi(lag_ind),'r')
legend('Rainfall Residuals','Typhi Residuals')
datetick('x','mmm-yy')
xlim([datenum(date_nepal(lag_ind(1),:))-7 datenum(date_nepal(lag_ind(end),:))+7])
subplot(1,2,2)
scatter(resid_rainfall,resid_typhi(lag_ind)) %change back to 1:716 for no lag
title('Rainfall vs Typhi Residuals')

%%
corr_resid=zeros(9,1);

for i=0:8
    corr_resid(i+1,1)=corr(resid_rainfall,resid_typhi(1+i:length(rainfall_nepalwk)+i,1));
    %(i+1,2)=corr(rainfall_nepalwk,paratyphi_nepal(1+i:length(rainfall_nepalwk)+i,1));
end

%%
year=zeros(length(date_nepal),1);
week=zeros(length(date_nepal),1);
week(1)=15;
for i=2:length(date_nepal)
    if date_nepal(i,1)==date_nepal(i-1,1)
        week(i,1)=week(i-1,1)+1;
    else
        week(i,1)=1;
        year(i,1)=1;
    end
end  
yri=find(year);

typhi_wkyr=zeros(52,14);
for j=2:13
    typhi_wkyr(:,j)=typhi_nepal(yri(j):yri(j)+51,1);
end
typhi_wkyr(15:52,1)=typhi_nepal(1:38,1);
typhi_wkyr(1:26,14)=typhi_nepal(yri(end):yri(end)+25,1);

%% -- check, maybe this is not necessary?
% rainday=(datenum([1997 1 1 0 0 0]):datenum([2010 12 28 0 0 0]))';
% date_rain=[1996 12 29 0 0 0];
% while date_rain(end,1)<2011
%     date_rain=[date_rain; datevec(datenum(date_rain(end,:))+7)];
% end
% 
% %figure out what's going on with rain dates here
% %get weekly rain totals
% rainfall_nepal=zeros(length(date_rain)-1,1);
% j=1;
% for i=1:length(rainday)
%     if rainday(i)<datenum(date_rain(j+1,:))
%         rainfall_nepal(j,1)=rainfall_nepal(j,1)+rainfall_nepalday(i,1);
%     else
%         j=j+1;
%         rainfall_nepal(j,1)=rainfall_nepalday(i,1);
%     end
% end
% 
% date_rain(end,:)=[];
%%
rainpar_nepal=fminsearch(@(p) seasfit(p,rainbywk_nepal,52),[50 40 .6])

figure
plot([rainbywk_nepal rainpar_nepal(1)+rainpar_nepal(2)*cos(2*pi*((1:52)'-rainpar_nepal(3)*52)/52)])

%% Do analysis for paratyphoid and rainfall

logparatyphi = log(paratyphi_nepal+1);
%% %% Fit harmonic regression

time=(15:length(date_nepal)+14)'; % data starts April 13 = week 15

mdl_paratyphi1=GeneralizedLinearModel.fit([time cos(2*pi*time/52.18) sin(2*pi*time/52.18)],logparatyphi,'linear','Distribution','poisson','Link','log');
mdl_paratyphi2=GeneralizedLinearModel.fit([time cos(2*pi*time/52.18) sin(2*pi*time/52.18) cos(2*pi*time/26.09) sin(2*pi*time/26.09)],logparatyphi,'linear','Distribution','poisson','Link','log');
%mdl_typhi3=GeneralizedLinearModel.fit([time cos(2*pi*time/52.18) sin(2*pi*time/52.18) cos(2*pi*time/26.09) sin(2*pi*time/26.09) cos(2*pi*time/13) sin(2*pi*time/13)],typhi_nepal,'linear','Distribution','poisson','Link','log');
%mdl_typhi4=GeneralizedLinearModel.fit([time cos(2*pi*time/156.5) sin(2*pi*time/156.5) cos(2*pi*time/52.18) sin(2*pi*time/52.18) cos(2*pi*time/26.09) sin(2*pi*time/26.09) cos(2*pi*time/13) sin(2*pi*time/13)],typhi_nepal,'linear','Distribution','poisson','Link','log');
%mdl_typhi5=GeneralizedLinearModel.fit([time cos(2*pi*time/156.5) sin(2*pi*time/156.5) cos(2*pi*time/52.18) sin(2*pi*time/52.18) cos(2*pi*time/26.09) sin(2*pi*time/26.09)],typhi_nepal,'linear','Distribution','poisson','Link','log');
%mdl_typhi6=GeneralizedLinearModel.fit([time cos(2*pi*time/521.8) sin(2*pi*time/521.8) cos(2*pi*time/156.5) sin(2*pi*time/156.5) cos(2*pi*time/52.18) sin(2*pi*time/52.18) cos(2*pi*time/26.09) sin(2*pi*time/26.09)],typhi_nepal,'linear','Distribution','poisson','Link','log');
mdl_paratyphi7=GeneralizedLinearModel.fit([time cos(2*pi*time/626.1) sin(2*pi*time/626.1) cos(2*pi*time/156.5) sin(2*pi*time/156.5) cos(2*pi*time/52.18) sin(2*pi*time/52.18) cos(2*pi*time/26.09) sin(2*pi*time/26.09)],logparatyphi,'linear','Distribution','poisson','Link','log');
%mdl_typhi8=GeneralizedLinearModel.fit([time cos(2*pi*time/(52.18*14)) sin(2*pi*time/(52.18*14)) cos(2*pi*time/156.5) sin(2*pi*time/156.5) cos(2*pi*time/52.18) sin(2*pi*time/52.18) cos(2*pi*time/26.09) sin(2*pi*time/26.09)],typhi_nepal,'linear','Distribution','poisson','Link','log');

mdl_paratyphi3=GeneralizedLinearModel.fit([time cos(2*pi*time/156.5) sin(2*pi*time/156.5) cos(2*pi*time/52.18) sin(2*pi*time/52.18)],logparatyphi,'linear','Distribution','poisson','Link','log');
mdl_paratyphi4=GeneralizedLinearModel.fit([time cos(2*pi*time/626.1) sin(2*pi*time/626.1) cos(2*pi*time/52.18) sin(2*pi*time/52.18)],logparatyphi,'linear','Distribution','poisson','Link','log');
mdl_paratyphi5=GeneralizedLinearModel.fit([time cos(2*pi*time/626.1) sin(2*pi*time/626.1) cos(2*pi*time/156.5) sin(2*pi*time/156.5) cos(2*pi*time/52.18) sin(2*pi*time/52.18)],logparatyphi,'linear','Distribution','poisson','Link','log');

%% display AICs
AICs_paratyphoid = [];
for i = [1:5,7]
   eval(['tmpMdl = mdl_paratyphi',num2str(i),';'])
   %disp(i);
   %disp(tmpMdl.ModelCriterion.AIC);
   AICs_paratyphoid = [AICs_paratyphoid;tmpMdl.ModelCriterion.AIC];
end

%%
%model 4 has lowest AIC again
%Y_typhi6=predict(mdl_typhi6);
%Y_typhi7=predict(mdl_typhi7);
%Y_typhi8=predict(mdl_typhi8);
Y_paratyphi4=max(0,exp(predict(mdl_paratyphi4))-1);

figure; hold on
plot(datenum(date_nepal),paratyphi_nepal,'b')
plot(datenum(date_nepal),Y_paratyphi4,'r')
%plot(datenum(date_nepal),Y_typhi6,'r')
%plot(datenum(date_nepal),Y_typhi7,'g')
%plot(datenum(date_nepal),Y_typhi8,'c')
datetick('x','mmm-yy')
xlim([datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7])
set(gca,'YLim',[0 80],'YTick',0:20:80)
ylabel({'No. of Paratyphi cases';'(per week)'})
legend('Paratyphoid Cases','Model Fit 4 (12, 1 year periods)'); %,'Model Fit 7 (12,3,1, .5 yr periods)')

%%
%best model is 4 again, take residuals and plot vs rainfall

resid_paratyphi=paratyphi_nepal./Y_paratyphi4;

figure
[ax,y1,y2]=plotyy(datenum(date_nepal),resid_paratyphi,datenum(date_rainwk),rainfall_nepalwk);
datetick('x','mmm-yy')
set(y1,'Color','r')
set(y2,'Color','b')
set(ax(1),'YColor','k','XLim',[datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7],'YTick',0:2:10)
set(ax(2),'YColor','k','XLim',[datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7],'XTick',[],'YTick',0:100:500)
ylabel('Paratyphoid Residuals (per week)')
set(get(ax(2),'YLabel'),'String','Rainfall (mm/wk)')
legend('Paratyphoid Residuals','Rainfall')

%% compare residuals

% define lag index vector
lag_ind= [1:716];
date_num_nepal = datenum(date_nepal);
figure; 
subplot(1,2,1); hold on;
plot(datenum(date_rainwk),resid_rainfall,'b')
%plot(datenum(date_nepal),resid_paratyphi,'r')
plot(date_num_nepal(lag_ind),resid_paratyphi(lag_ind),'r')
legend('Rainfall Residuals','Paratyphi Residuals')
datetick('x','mmm-yy')
xlim([datenum(date_nepal(lag_ind(1),:))-7 datenum(date_nepal(lag_ind(end),:))+7])
subplot(1,2,2)
scatter(resid_rainfall,resid_paratyphi(lag_ind))
title('Rainfall vs Paratyphi Residuals')

%% Calculate correlation between rainfall and typhi/paratyphi cases and different lags (0-8 weeks)

corr_rainTPresid=zeros(9,1);
corr_rainALLresid=zeros(9,1);

for i=0:8
    corr_rainTPresid(i+1,1)=corr(rainfall_nepalwk,resid_typhi(1+i:length(rainfall_nepalwk)+i,1));
    corr_rainTPresid(i+1,2)=corr(rainfall_nepalwk,resid_paratyphi(1+i:length(rainfall_nepalwk)+i,1));
    corr_rainALLresid(i+1,1)=corr(resid_rainfall,resid_typhi(1+i:length(rainfall_nepalwk)+i,1));
    corr_rainALLresid(i+1,2)=corr(resid_rainfall,resid_paratyphi(1+i:length(rainfall_nepalwk)+i,1));

end

%% compare residuals-- not with rainfall, just to test

figure; 
subplot(1,2,1); hold on;
plot(datenum(date_nepal),resid_typhi,'b')
plot(datenum(date_nepal),resid_paratyphi,'r')
legend('Typhoid Residuals','Parayphi Residuals')
datetick('x','mmm-yy')
xlim([datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7])
subplot(1,2,2)
scatter(resid_typhi,resid_paratyphi)
title('Typhoid vs Parayphi Residuals')

%% Lag Analysis Figures
Y_rainfall6=predict(mdl_rainfall6);
    resid_rainfall=rainfall_nepalwk./Y_rainfall6;
    
for i = 0:8
    % define lag index vector
    lag_ind= [1+i:716+i];
    date_num_nepal = datenum(date_nepal);
    figure; 
    subplot(1,2,1); hold on;
    plot(datenum(date_rainwk),resid_rainfall,'b')
    %plot(datenum(date_nepal),resid_typhi,'r')
    plot(date_num_nepal(lag_ind),resid_typhi(lag_ind),'r')
    legend('Rainfall Residuals','Typhi Residuals')
    datetick('x','mmm-yy')
    xlim([datenum(date_nepal(lag_ind(1),:))-7 datenum(date_nepal(lag_ind(end),:))+7])
    subplot(1,2,2)
    scatter(resid_rainfall,resid_typhi(lag_ind)) %change back to 1:716 for no lag
    title(['Rainfall vs Typhi Residuals, ',num2str(i),' week lag'])


    figure; 
    subplot(1,2,1); hold on;
    plot(datenum(date_rainwk),resid_rainfall,'b')
    %plot(datenum(date_nepal),resid_paratyphi,'r')
    plot(date_num_nepal(lag_ind),resid_paratyphi(lag_ind),'r')
    legend('Rainfall Residuals','Paratyphi Residuals')
    datetick('x','mmm-yy')
    xlim([datenum(date_nepal(lag_ind(1),:))-7 datenum(date_nepal(lag_ind(end),:))+7])
    subplot(1,2,2)
    scatter(resid_rainfall,resid_paratyphi(lag_ind))
    title(['Rainfall vs Paratyphi Residuals, ',num2str(i),' week lag'] )
end