%% 1. Arrange the data in a 26x21 biweek by year matrix, calculate the mean 
% and range for each year, then scatter plot the range vs mean

figure
scatter(mean(typhi_wkyr),range(typhi_wkyr))
xlabel('Mean')
ylabel('Range')

[r,pval]=corr(mean(typhi_wkyr)',range(typhi_wkyr)')

%% 2. Transform the data by taking the log of the number of cases, then 
% repeat question 1 above. Plot the log-transformed data alongside the
% mean-range plot.  Does it appear as though the requirements for
% stationarity are met?

logcases=log(typhi_nepal+1);
logcases_byyr=log(typhi_wkyr+1); 

figure
subplot(1,2,1)
plot(datenum(date_nepal),logcases) 
datetick('x','mmm-yy')
xlim([datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7])
ylabel('Log(Typhoid cases)')
subplot(1,2,2)
scatter(mean(logcases_byyr)',range(logcases_byyr)')
xlabel('Mean')
ylabel('Range')

[r,pval]=corr(mean(logcases_byyr)',range(logcases_byyr)')


%% 3. Calculate and plot the autocorrelation and partial autocorrelation functions.
figure
subplot(2,1,1)
autocorr(logcases,52*6,0,4) % Calculates ACF up to 6 years, with a significance cutoff line plotted at 4 standard deviations (~90% CI)
title('ACF')
subplot(2,1,2)
parcorr(logcases,52*6,0,4) % Calculates PACF up to 6 years, with a significance cutoff line plotted at 4 standard deviations (~90% CI)
title('PCF')

%% 4. Define and fit and ARIMA model 
model=arima('Constant',0,'D',0,'Seasonality',52,'ARLags',1,'MALags',1,'SARLags',[52 104 156])
%model=arima('Constant',0,'D',0,'Seasonality',0,'ARLags',[1:3 52],'MALags',1:2)%,'SARLags',52)%:52:260)
[fit,~,LL]=estimate(model,logcases)

%% 5. Calculate and plot the model residuals, along with their ACF and PACF

resid_arima=infer(fit,logcases);

figure
subplot(3,1,1)
plot(datenum(date_nepal),resid_arima)
datetick('x','mmm-yy')
xlim([datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7])
ylabel('Model residuals')
subplot(3,1,2)
autocorr(resid_arima,52*6,0,4)
title('Residuals ACF')
subplot(3,1,3)
parcorr(resid_arima,52*6,0,4)
title('Residuals PACF')

%% 6. Fit the ARIMA model to the first 10+ years of data, then 
% use the fitted model to predict the last 2.5 years of the data.

Y1 = logcases(1:560);
fit1 = estimate(model,Y1);
Yfit1 = logcases(1:560)-infer(fit1,Y1);
Ypredict1 = forecast(fit1,length(logcases)-560,'Y0',Y1);

figure
subplot(1,2,1)
hold on
plot(datenum(date_nepal),logcases,'k')
plot(datenum(date_nepal(1:560,:)),Yfit1,'r')
plot(datenum(date_nepal(561:end,:)),Ypredict1,'r')
plot([datenum(date_nepal(560,:)) datenum(date_nepal(560,:))],[0 4.5],'--r')
datetick('x','mmm-yy')
xlim([datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7])
ylabel('log(cases)')
legend('Observed','Forecast','Location','NorthWest')
subplot(1,2,2)
hold on
plot(datenum(date_nepal),exp(logcases)-1,'k')
plot(datenum(date_nepal(1:560,:)),exp(Yfit1)-1,'r')
plot(datenum(date_nepal(561:end,:)),exp(Ypredict1)-1,'r')
plot([datenum(date_nepal(560,:)) datenum(date_nepal(560,:))],[0 80],'--r')
datetick('x','mmm-yy')
xlim([datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7])

%%
figure
[ax,y1,y2]=plotyy(datenum(date_nepal),resid_arima,datenum(date_rainwk),rainfall_nepalwk);
datetick('x','mmm-yy')
set(y1,'Color','r')
set(y2,'Color','b')
set(ax(1),'YColor','k','XLim',[datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7],'YTick',0:20:100)
set(ax(2),'YColor','k','XLim',[datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7],'XTick',[],'YTick',0:100:500)
ylabel('Typhoid cases (per week)')
set(get(ax(2),'YLabel'),'String','Rainfall (mm/wk)')
legend('ARIMA Residuals','Rainfall')

corr_rainarimaresid=zeros(9,1);
for i=0:8
    corr_rainarimaresid(i+1,1)=corr(rainfall_nepalwk,resid_arima(1+i:length(rainfall_nepalwk)+i,1));
    %(i+1,2)=corr(rainfall_nepalwk,paratyphi_nepal(1+i:length(rainfall_nepalwk)+i,1));
end


%% 3. Calculate and plot the autocorrelation and partial autocorrelation functions for log rainfall.
figure
subplot(2,1,1)
autocorr(lograinfall,52*6,0,4) % Calculates ACF up to 6 years, with a significance cutoff line plotted at 4 standard deviations (~90% CI)
title('ACF')
subplot(2,1,2)
parcorr(lograinfall,52*6,0,4) % Calculates PACF up to 6 years, with a significance cutoff line plotted at 4 standard deviations (~90% CI)
title('PCF')

%% 4. Define and fit and ARIMA model 
%model=arima('Constant',0,'D',0,'Seasonality',52,'ARLags',1:2,'SARLags',52:52:104)
model=arima('Constant',0,'D',0,'Seasonality',52,'ARLags',1,'MALags',1,'SARLags',[52 104 156])%:52:260)
[fit,~,LL]=estimate(model,lograinfall)

%% 5. Calculate and plot the model residuals, along with their ACF and PACF

resid_rainarima=infer(fit,lograinfall);

figure
subplot(3,1,1)
plot(datenum(date_rainwk),resid_rainarima)
datetick('x','mmm-yy')
xlim([datenum(date_rainwk(1,:))-7 datenum(date_rainwk(end,:))+7])
ylabel('Model residuals')
subplot(3,1,2)
autocorr(resid_rainarima,52*6,0,4)
title('Residuals ACF')
subplot(3,1,3)
parcorr(resid_rainarima,52*6,0,4)
title('Residuals PACF')

%% 6. Fit the ARIMA model to the first 10 years of rainfall notifications, then 
% use the fitted model to predict the last 2.5 years of the data.

R1 = lograinfall(1:560);
fit1r = estimate(model,R1);
Rfit1 = lograinfall(1:560)-infer(fit1r,R1);
Rpredict1 = forecast(fit1r,length(lograinfall)-560,'Y0',R1);

figure
subplot(1,2,1)
hold on
plot(datenum(date_rainwk),lograinfall,'b')
plot(datenum(date_rainwk(1:560,:)),Rfit1,'r')
plot(datenum(date_rainwk(561:end,:)),Rpredict1,'r')
plot([datenum(date_rainwk(560,:)) datenum(date_rainwk(560,:))],[0 4.5],'--r')
datetick('x','mmm-yy')
xlim([datenum(date_rainwk(1,:))-7 datenum(date_rainwk(end,:))+7])
ylabel('log(rainfall)')
legend('Observed','Forecast','Location','NorthWest')
subplot(1,2,2)
hold on
plot(datenum(date_rainwk),exp(lograinfall)-1,'k')
plot(datenum(date_rainwk(1:560,:)),exp(Rfit1)-1,'r')
plot(datenum(date_rainwk(561:end,:)),exp(Rpredict1)-1,'r')
plot([datenum(date_rainwk(560,:)) datenum(date_rainwk(560,:))],[0 80],'--r')
datetick('x','mmm-yy')
xlim([datenum(date_rainwk(1,:))-7 datenum(date_rainwk(end,:))+7])

%%
figure; 
subplot(1,2,1); hold on;
plot(datenum(date_rainwk),resid_rainarima,'b')
plot(datenum(date_nepal),resid_arima,'r')
datetick('x','mmm-yy')
xlim([datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7])
legend('Rainfall','Typhoid')
title('Residuals over Time')
subplot(1,2,2)
scatter(resid_rainarima,resid_arima(1:716),'ok')
title('Rainfall ARIMA Residuals vs Typhoid ARIMA Residuals')

corr_arimaresid=zeros(9,1);

for i=0:8
    corr_arimaresid(i+1,1)=corr(resid_rainarima,resid_arima(1+i:length(rainfall_nepalwk)+i,1));
end
