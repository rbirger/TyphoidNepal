clear B u mu omega r rC alpha delta theta xi gamma epsilon betaw betap q lag al v1 massvacc veff wv tcarr immig; 
global B u mu omega r rC alpha delta theta xi gamma epsilon betaw betap q lag al v1 massvacc veff wv tcarr san immig; 

%DEMOGRAPHIC PARAMETERS

%age = [2.5 7.5 15:10:85];
age = 2.5:5:82.5;
%agecat = {'0-4','5-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
agecat = {'0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79','80+'};
al = length(age); %1; %         %number of age groups
%agep = [260*ones(1,2) 520*ones(1,8)]; %1; %
agep = ones(1,al);
agep = agep/sum(agep); %proportion of population in each age group

%u = [1/260*ones(1,4) 1/520*ones(1,6) 0];  %0;  %rate of aging out of age group
u = (1/260)*ones(1,al);
mu = log(1.03)/52;   %natural death rate (per week)

t0=1040; %2600; %length of burn-in period
tspan=t0+1040; %length(C_nepal)+15; %

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
%r = 0.7;        %relative infectiousness of second infections
%alpha = .005;    %proportion of primary infections that die 
delta = .25;    %rate of recovery (1/weeks)
theta = .03;    %.03 proportion of infecteds who become carriers
xi = .33;         %rate of decay of infectious particles from water

%Estimated parameters

par1.san = ones(tspan,1);       %effect of sanitation 
par1.betap = .2*ones(al); %.0000000791*ones(al)*sum(N0); %person-to-person transmission rate
par1.betaw = .0000000025;  %.000000237;  %exp(-san*(1:tspan-t0-813)') transmission rate from water
par1.gamma = 1;  %rate of shedding into water
par1.q = 10;  %7.5 %amplitude of seasonal forcing
par1.lag = 0;       %timing of peak transmission
par1.r = 0.02;   %0.7;  %      %relative infectiousness of non-primary infections
par1.rC = 0.02;  %0.7;  %       %relative infectiousness of carriers
par1.epsilon = 0.022;  %0.02;   %fraction of recovereds who revert to fully susceptible state
par1.alpha = 0.5;    %mortality/reporting rate

par2.san = ones(tspan,1);       %effect of sanitation 
par2.betap = .00000008*ones(al)*sum(N0); %person-to-person transmission rate
par2.betaw = .000000237*.000129;  %exp(-san*(1:tspan-t0-813)') transmission rate from water
par2.gamma = 1;  %rate of shedding into water
par2.q = 121.17;  %amplitude of seasonal forcing
par2.lag = 0;       %timing of peak transmission
par2.r = 0.7;  %      %relative infectiousness of non-primary infections
par2.rC = 0.7;  %       %relative infectiousness of carriers
par2.epsilon = 0.022;   %fraction of recovereds who revert to fully susceptible state
par2.alpha = 0.5;    %mortality/reporting rate

par=par2; %pfit2;

san = par.san; %[ones(t0+tvacc,1); .05*ones(tspan-t0-tvacc,1)]; %
betap = par.betap;
betaw = par.betaw;
gamma = 1; %par.gamma;
q = par.q;
lag = 0; %par.lag;
r = par.r;
rC = par.rC;
epsilon = par.epsilon;
alpha = par.alpha;

betaw = betaw*(1+q*rainbywk_nepal);
%betaw = betaw*(1+q*[mean(rainfall_nepalwk)*ones(t0-52,1); rainbywk_nepal; rainbywk_nepal(1:15); rainfall_nepalwk; rainbywk_nepal(1:26)]);

R0=(1/(delta+mu))*(betap(1) + mean(betaw)*sum(N0)*gamma/xi)*(1 + rC*delta*theta/mu);
R0p=betap(1)/(delta+mu)*(1 + rC*delta*theta/mu);
R0w=(mean(betaw)*sum(N0)*gamma)/(xi*(delta+mu))*(1 + rC*delta*theta/mu);
R0c=(rC*delta*theta/mu)/(1+rC*delta*theta/mu);

%for sim=1:100
tcarr = zeros(tspan,1); %0.1*.25*[zeros(t0+260,1); ones(tspan-t0-260,1)]; %
%tcarr(t0+260,1)=0.5;
%Vaccination
for v=1 %[1 6] %1:11
vcov=(v-1)/10; %0.8; %vaccination coverage
veff=0.7; %vaccine efficacy
wv=1/156; %waning of vaccine immunity
mass=1;
agevacc=2; %age at vaccination (1=birth, 2=5yrs, 3=10yrs, etc)
tvacc=260; %week of initiating vaccination
v1=zeros(tspan,al);
massvacc=zeros(tspan,al);
if mass==1
    massvacc(t0+tvacc-4,:)=vcov*veff;
else
    massvacc(t0+tvacc-4,2:3)=vcov*veff;
    v1(:,agevacc)=[zeros(t0+tvacc,size(agevacc,2)); ones(tspan-t0-tvacc,size(agevacc,2))]*vcov*veff; 
end

%Initial conditions, S1, I1, R, S2, I2, C, W 
y0 = [round(.88*N0')-10;    %S1
    10*ones(al,1);                  %I1
    round(0*N0');         %R
    round(.1*N0')-10;      %S2
    10*ones(al,1);                  %I2
    round(.02*N0');         %C
    zeros(al,1);    %W 
    zeros(al,1);
    10];            %V

%Differential equations
%options=odeset('NonNegative',1:length(y0));
options = odeset('RelTol',1e-6) ;
[t,y] = ode45(@typhoid_ode_nepalv,1:tspan,y0,options);

%Discard burn-in period
t(1:t0)=[]; 
y(1:t0,:)=[];

%%
lambdap=((y(:,al+1:2*al)+r*y(:,4*al+1:5*al)+rC*y(:,5*al+1:6*al))./(sum(y(:,1:8*al),2)*ones(1,al)))*betap; %
lambdaw=betaw(rem(round(t),52)+1,1).*y(:,8*al+1); %./sum(y(:,1:6*al),2)
W=y(:,8*al+1);
pop=y(:,1:al);
for i=2:8
    pop=pop+y(:,(i-1)*al+1:i*al);
end
C=y(:,5*al+1:6*al);

Dv=delta*alpha*y(:,al+1:2*al); % 
%D=poissrnd(D);
agedistm=Dv./(sum(Dv,2)*ones(1,al)); %D(417:677,:)./(sum(D(417:677,:),2)*ones(1,al)); %
A=agedistm*age';

%cumulative proportion of individuals vaccinated in each age group
V=y(:,6*al+1:7*al)./pop;
doses=(sum((v1(t0+1:tspan,:).*pop.*(ones(tspan-t0,1)*u)),2)./sum(pop,2)+massvacc(t0+1:tspan,1))*100000/veff;
if vcov==0
    Dnovaccv=Dv;
end
directeffv(:,v)=sum(V.*Dnovaccv,2)*veff; %*(1-epsilon)

%preduct_rand1(v,sim)=1-sum(sum(D(tvacc+1:tvacc+52,:),2))/sum(sum(Dnovacc(tvacc+1:tvacc+52,:),2));
preductv_y1(v,1)=1-sum(sum(Dv(tvacc+1:tvacc+52,:),2))/sum(sum(Dnovaccv(tvacc+1:tvacc+52,:),2));
pdirectv_y1(v,1)=sum(directeffv(tvacc+1:tvacc+52,v))/sum(sum(Dnovaccv(tvacc+1:tvacc+52,:),2));
pcoveffv_y1(v,1)=vcov*veff;

end
%end
%%
figure
subplot(2,1,1)
hold on
plot(1/52:1/52:20,sum(Dnovaccv,2),':b','LineWidth',2)
plot(1/52:1/52:20,sum(Dnovaccv,2)-directeffv(:,v),'g','LineWidth',2)
plot(1/52:1/52:20,sum(Dv,2),'r','LineWidth',2)
plot(1/52:1/52:20,sum(Dnovaccv,2),':b','LineWidth',2)
%ylim([0 50])
%subplot(3,1,2)
%hold on
%plot(sum(y(:,1:6*al),2))
%plot(sum(N,2),'r')
subplot(2,1,2)
bar(mean(agedistm)') 
set(gca,'XTick',1:2:al,'XTickLabel',agecat(1:2:al))
%plot(1/52:1/52:20,sum(pop,2))

%%
%figure
%hold on
%plot(15:756,C_nepal(:,1))
%plot(15:756,Dnepalfit1,'r')
%plot(15:756,Dnepalfit2,'g')

