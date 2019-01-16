%% Fourier analysis

tmax=length(paratyphi_nepal); %length of the time series
dt=1/52; %time step (in years)

ym=abs(fft(paratyphi_nepal)); %absolute value of the fast Fourier transform for the time series 

f=(0:tmax-1)/tmax; %frequency of oscillations (in bi-weeks)
T=dt./f; %period of oscillations (in years) = 1/frequency

figure
bar(T(2:round(tmax/2)+1),ym(2:round(tmax/2)+1,:),'EdgeColor','k')
xlim([0 10])
xlabel('Period (years)','FontSize',12)
ylabel('Strength','FontSize',12)

%% Wavelet analysis  

% First, log-transform the data and scale it to have a mean=0 and variance=1 
logtyphi=log(paratyphi_nepal+1); 
lograinfall=log(rainfall_nepalwk+1);

tsdata=(logtyphi - mean(logtyphi))/std(logtyphi); %lograinfall; %

n = length(tsdata);  % number of time points
dt = 1/52;   % frequency of observations
pad = 1;      % pad the time series with zeroes (recommended)
dj = 0.25;    % spacing between scales (this will do 4 sub-octaves per octave)
s0 = 13*dt;    % this says start at a scale of 13 weeks
j1 = 6/dj;    % this says do 6 powers-of-two with dj sub-octaves each (T=2^0 to 2^5=32)
lag1 = corr(tsdata(1:n-1),tsdata(2:n));  % lag-1 autocorrelation for red noise background
mother = 'Morlet';

% Wavelet transform:
[wave,period,scale,coi] = wavelet(tsdata,dt,pad,dj,s0,j1,mother);
power = (abs(wave)).^2 ;        % compute wavelet power spectrum
power = power*dt;
% Significance levels: (variance=1 for the normalized data)
[signif,fft_theor] = wave_signif(1.0,dt,scale,0,lag1); %,-1,-1,mother);
sig95 = (signif')*(ones(1,n));  % expand signif --> (J+1)x(N) array
sig95 = power ./ sig95;         % where ratio > 1, power is significant

% Global wavelet spectrum & significance levels:
global_ws = (std(tsdata)^2)*(sum(power')/n);   % time-average over all times
dof = n - scale;  % the -scale corrects for padding at edges
global_signif = wave_signif(std(tsdata)^2,dt,scale,1,lag1,-1,dof,mother);

%create time vector
time = dt*[1:tmax]';

figure
%--- Plot time series
subplot(3,4,1:3)
plot(time(1:n),tsdata)
xlim([time(1,:) time(n,:)])
xlabel('Time (years)')
ylabel({'Normalized log number of'; 'cases (per week)'},'FontSize',10)

%--- Contour plot (filled) of wavelet power spectrum
subplot(3,4,5:7)
levels = -8:10; %specifies log2 of where contour levels should be
contourf(time(1:n),log2(period),max(log2(power),-8),-8:10,'LineColor','none');  %log2(levels) *** use 'contourf' for "filled" contour plot
%imagesc(time,log2(period),log2(power));  %*** uncomment for 'image' plot
% 95% significance contour, levels at -99 (fake) and 1 (95% signif)
hold on
contour(time(1:n),log2(period),abs(sig95),[-99,1],'k','LineWidth',2);
plot(time(1:n),log2(coi),'w','LineWidth',2) % cone-of-influence, anything "above" is dubious
%xlabel('Year')
ylabel('Period (years)','FontSize',10)
title('Wavelet Power Spectrum','FontSize',10)
datetick('x','yyyy')
set(gca,'FontSize',10)
set(gca,'XLim',[time(1,:) time(n,:)])
set(gca,'YLim',[-2 4],'YTick',-2:4,'YTickLabel',[0.25 0.5 1 2 4 8 16])

%--- Plot global wavelet spectrum
subplot(3,4,8)
plot(global_ws,log2(period))
hold on
plot(abs(global_signif),log2(period),':k')
xlabel('Power','FontSize',10)
ylabel('Period (years)','FontSize',10)
title('Global Wavelet Spectrum','FontSize',10)
set(gca,'YLim',[-1 4],'YTick',-1:4,'YTickLabel',[0.5 1 2 4 8 16])
set(gca,'XLim',[0 5])

%--- Plot dominant epidemic period and compare to the birth rate
[~,xm]=max(power);
subplot(3,4,9:11)
plot(time(1:n),period(xm));
set(gca,'XLim',[time(1,:) time(n,:)])
ylabel('Dominant period')

%%
% 7. Calculate the phase angle for each scale and point in the time series,
% then plot the phase angle (in degrees) corresponding to the 1-year 
% period, as in Fig. 2d of Grenfell et al (2001).


phase=angle(wave);

figure
plot(time(1:n),phase(9,:)*180/pi)
xlim([time(1,:) time(n,:)])
xlabel('Time (years)')
ylabel('Phase angle (degrees)')
