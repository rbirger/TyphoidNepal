clear B u mu omega r alpha delta theta xi gamma epsilon betaw betap q lag al v1 massvacc immig san t1 t2 rt; % noise; rC 
global B u mu omega r alpha delta theta xi gamma epsilon betaw betap q lag al v1 massvacc immig san t1 t2 rt; % noise; rC 

%DEMOGRAPHIC PARAMETERS

%age = [2.5 7.5 15:10:85];
age = [2.5:5:82.5];
%agecat = {'0-4','5-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
agecat = {'0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79','80+'};
al = length(age); %1; %         %number of age groups
%agep = [260*ones(1,2) 520*ones(1,8)]; %1; %
agep = ones(1,al);
agep = agep/sum(agep); %proportion of population in each age group

%u = [1/260*ones(1,4) 1/520*ones(1,6) 0];  %0;  %rate of aging out of age group
u = (1/260)*ones(1,al);
mu = log(1.03)/52;   %natural death rate (per week)

t0=round(52.18*47); %2600; %length of burn-in period
tspan=t0+length(typhi_nepal)+15; %1040; %

datesim=(datenum([1950 1 1 0 0 0]):7:datenum([1950 1 1 0 0 0])+7*(tspan-1))';

N = 100000*agep;    %93859 Population size 
N0=N(1,:); %/exp((B(1)/N(1)-mu)*t0);          
brate = .03;  %Crude birth rate (per year)
B = zeros(size(N));
B(:,1) = log(1+brate)/52*sum(N0,2); %*ones(length(B),1);    %Number of births per week

immig = zeros(tspan,al);
%immig(t0+156:t0+208,4:7)=50;


%MODEL PARAMETERS

%Fixed parameters

omega = 1/104;  %rate of waning immunity (1/weeks)
alpha = .001;    %proportion of primary infections that die 
delta = .25;    %rate of recovery (1/weeks)
theta = [.003*ones(4,1); .021*ones(2,1); .044*ones(2,1); .088*ones(2,1); .101*ones(2,1); .078*ones(5,1)]; %.03;    %proportion of infecteds who become carriers
xi = .33;         %rate of decay of infectious particles from water

gamma = 1; %
san = ones(tspan,1); %
epsilon = 0; 


%Estimated parameters

par0=[1; 1; .4; .2; .01; 200; 300; 1.3]; %[1.5; 0.2; 0; .05; .01; 52; 104; 1.3];
par1=[.8; 2; .3; .2; .1; 200; 300; 3];
par=par1; %par1_nepal; %

R0p = par(1); %.5*par(1); %
R0w = par(2); %.5*par(1); % 
q = par(3); %rainpar(2)/rainpar(1); %
lag = rainpar_nepal(3)*52.18; %0; %       %timing of peak transmission
rep = par(4); 
r = par(5); %.1; % 
rC = r; 

t1 = par(6)+t0; %t0; %
t2 = par(7)+t0; %t0; %
rt = par(8); %1; %

ximm = 1000; %par(6);
wx = 1/52; %1/par(7); %rate of waning cross-immunity from NTS infection


betap = R0p*(delta+mean(mu))/(1 + rC*delta*(agep*theta)/mean(mu))*ones(al); %person-to-person transmission rate
%betaw = R0w/(N1*gamma)*(xi*(delta+mean(mu)))/(1 + rC*delta*(agep*theta)/mean(mu)); %/(q*mean(rainfall)); %transmission rate from water (density dependent)
betaw = (R0w/gamma)*(xi*(delta+mean(mu)))/(1 + rC*delta*(agep*theta)/mean(mu)); %/(q*mean(rainfall)); %transmission rate from water (frequency dependent)

betaw=betaw*ones(tspan,1);
%betaw = betaw*(1+q*rainbywk_nepal);
%betaw = betaw*(1+q*[mean(rainfall_nepalwk)*ones(t0-52,1); rainbywk_nepal; rainbywk_nepal(1:15); rainfall_nepalwk; rainbywk_nepal(1:26)]);

trans=ones(tspan,1);
dur=1/delta*ones(tspan,1);
for i=1:tspan
    if i>=t1 && i<t2
        trans(i)=1+(rt-1)*(i-t1)/(t2-t1);
        %dur(i,1)=dur(i,1)*(1+(rt-1)*(i-t1)/(t2-t1));
    elseif i>=t2
        trans(i)=rt;
        %dur(i,1)=dur(i,1)*rt;
    end
end

R0=(1./(1./dur+mean(mu))).*(betap(1) + mean(betaw)*gamma/xi).*trans.*(1 + rC*(1./dur)*(agep*theta)/mean(mu)); %frequency-dependent transmission
%R0=(1/(delta+mu))*(betap(1) + mean(betaw)*sum(N0)*gamma/xi)*(1 + rC*delta*(agep*theta)/mu)*trans;
R0c=rC*delta*(agep*theta)/mu;

%for sim=1:100
%Vaccination
for v=1 %1:11 %
vcov=(v-1)/10; %vaccination coverage
veff=0.7; %vaccine efficacy
mass=0;
agevacc=2; %age at vaccination (1=birth, 2=5yrs, 3=10yrs, etc)
tvacc=260; %week of initiating vaccination
v1=zeros(tspan,al);
massvacc=zeros(tspan,1);
if mass==1
    massvacc(t0+tvacc-4,1)=vcov*veff;
else
    v1(:,agevacc)=[zeros(t0+tvacc,size(agevacc,2)); ones(tspan-t0-tvacc,size(agevacc,2))]*vcov*veff; 
end

%Initial conditions, S1, I1, R, S2, I2, C, W 
y0 = [round(.88*N0')-10;    %S1
    10*ones(al,1);                  %I1
    round(0*N0');         %R
    round(.1*N0')-10;      %S2
    10*ones(al,1);                  %I2
    round(.02*N0');         %C
    10;            %W 
    zeros(al,1)];            %V

%Differential equations
%options=odeset('NonNegative',1:length(y0));
options = odeset('RelTol',1e-6) ;
[t,y] = ode45(@typhoid_ode_nepal,1:tspan,y0,options);

%Discard burn-in period
t(1:t0)=[]; 
y(1:t0,:)=[];

%%
lambdap=((y(:,al+1:2*al)+r*y(:,4*al+1:5*al)+rC*y(:,5*al+1:6*al))./(sum(y(:,1:6*al),2)*ones(1,al)))*betap; %
lambdaw=betaw(rem(round(t),52)+1,1).*y(:,6*al+1); %./sum(y(:,1:6*al),2)

C=y(:,5*al+1:6*al);
W=y(:,6*al+1);
X=y(:,6*al+2:7*al+1);

D=rep*delta*y(:,al+1:2*al);  
%D=poissrnd(D);

agedistm=D./(sum(D,2)*ones(1,al)); %D(417:677,:)./(sum(D(417:677,:),2)*ones(1,al)); %
A=agedistm*age';
pop=y(:,1:al);
for i=2:6
    pop=pop+y(:,(i-1)*al+1:i*al);
end

%cumulative proportion of individuals vaccinated in each age group
%V=y(:,6*al+2:7*al+1)./pop;
%doses=(sum((v1(t0+1:tspan,:).*pop.*(ones(tspan-t0,1)*u)),2)./sum(pop,2)+massvacc(t0+1:tspan,1))*100000/veff;
%if vcov==0
%    Dnovacc=D;
%end
%directeff(:,v)=sum(V.*Dnovacc,2)*veff; %*(1-epsilon)

%preduct_rand1(v,sim)=1-sum(sum(D(tvacc+1:tvacc+52,:),2))/sum(sum(Dnovacc(tvacc+1:tvacc+52,:),2));
%preduct_y1(v,2)=1-sum(sum(D(tvacc+1:tvacc+52,:),2))/sum(sum(Dnovacc(tvacc+1:tvacc+52,:),2));
%pdirect_y1(v,2)=sum(directeff(tvacc+1:tvacc+52,v))/sum(sum(Dnovacc(tvacc+1:tvacc+52,:),2));
%pcoveff_y1(v,2)=vcov*veff;

end
%end
%%
figure
subplot(2,1,1)
hold on
plot(datenum(date_nepal),typhi_nepal,'r')
[ax,y1,y2]=plotyy(datesim(t0+1:end),sum(D,2),datesim(t0+1:end),R0(t0+1:end,1));
set(ax(1),'XLim',[datesim(t0) datenum(datesim(end))+7],'XTick',[],'YColor','k')
%set(ax(1),'XLim',[datenum([1995 12 31 0 0 0]) datenum([2016 1 1 0 0 0])],'XTick',[],'YColor','k')
set(ax(2),'XLim',[datesim(t0) datenum(datesim(end))+7],'XTick',[],'YColor','k')
datetick('x','mmm-yy')
xlim([datenum(date_nepal(1,:))-7 datenum(date_nepal(end,:))+7])
%plot(1/52:1/52:20,sum(Dnovacc,2),':b','LineWidth',2)
%plot(1/52:1/52:20,sum(Dnovacc,2)-directeff(:,v),'g','LineWidth',2)
%plot(1/52:1/52:20,sum(D,2),'r','LineWidth',2)
%plot(1/52:1/52:20,sum(Dnovacc,2),':b','LineWidth',2)
%ylim([0 50])
subplot(2,1,2)
bar(mean(agedistm)') 
set(gca,'XTick',1:2:al,'XTickLabel',agecat(1:2:al))
%plot(1/52:1/52:20,sum(pop,2))

%%

