function [LL,D,R0]=typhoidfit_nepal(unkp,data,Npop,Brate,vargin)

global B u mu omega r alpha delta theta xi gamma epsilon betaw betap q lag al v1 massvacc immig san t1 t2 rt; % noise; rC 

if ~isempty(vargin)
    rain=vargin;
    %r=vargin;
end

priorLL=0;
priorpars=[0 10; 0 10; 0 1; 0 1; 0 1; 0 length(data); 0 length(data); 1 10];
%priorpars=[0 10; 0 10; 0 1; 0 1; 0 length(data); 0 length(data); 1 10];
%priorpars=[0 10; 0 10; 0 1; 0 1; 0 1; 1 sum(Npop); 1 1000];
%priorpars=[0 10; 0 10; 0 1; 0 1; 0 1];
for i=1:length(unkp)
    priorLL=priorLL+log(unifpdf(unkp(i),priorpars(i,1),priorpars(i,2)));
end
if priorLL>-100000

%DEMOGRAPHIC PARAMETERS

age = 2.5:5:82.5;
al = length(age); %1; %         %number of age groups
agep = ones(1,al);
agep = agep/sum(agep); %proportion of population in each age group

u = (1/260)*ones(1,al);
mu = log(1.03)/52;   %natural death rate (per week)

t0=round(52.18*47)+15; %2600; %length of burn-in period
tspan=t0+length(data); %1040; %

N = Npop*agep;    %93859 Population size 
N0=N(1,:); %/exp((B(1)/N(1)-mu)*t0);          
B = zeros(size(N));
B(:,1) = log(1+Brate)/52*sum(N0,2); %*ones(length(B),1);    %Number of births per week

immig = zeros(tspan,al);


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

par=unkp;

R0p = par(1); 
R0w = par(2); 
q = par(3); 
lag = 28.88;
rep = par(4); 
r = par(5); 
rC = r; 

t1 = par(6)+t0; %t0; %
t2 = par(7)+t0; %t0; %
rt = par(8); %1; %

%timm0 = 156;
%timmf = 208;
%imm = 0;

%immig(t0+round(timm0):t0+round(timmf),1) = imm;


betap = R0p*(delta+mean(mu))/(1 + rC*delta*(agep*theta)/mean(mu))*ones(al); %person-to-person transmission rate
%betaw = R0w/(N1*gamma)*(xi*(delta+mean(mu)))/(1 + rC*delta*(agep*theta)/mean(mu)); %/(q*mean(rainfall)); %transmission rate from water (density dependent)
betaw = (R0w/gamma)*(xi*(delta+mean(mu)))/(1 + rC*delta*(agep*theta)/mean(mu)); %/(q*mean(rainfall)); %transmission rate from water (frequency dependent)

%yrs=floor(length(rain)/52);
%rainbywk=rain(1:52);
%for i=1:yrs-1
%    rainbywk=rainbywk+rain(52*i+1:52*(i+1));
%end
%rainbywk=rainbywk/yrs;

betaw=betaw*ones(tspan,1);
%betaw=betaw*(1+q*[mean(rain)*ones(t0-52,1); rainbywk; rain; rainbywk]);

trans=ones(tspan,1);
%dur=1/delta*ones(tspan,1);
for i=1:tspan
    if i>=t1 && i<t2
        trans(i)=1+(rt-1)*(i-t1)/(t2-t1);
        %dur(i,1)=dur(i,1)*(1+(rt-1)*(i-t1)/(t2-t1));
    elseif i>=t2
        trans(i)=rt;
        %dur(i,1)=dur(i,1)*rt;
    end
end
R0=(par(1)+par(2))*trans;

v1=zeros(tspan,al);
massvacc=zeros(tspan,1);

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

options = odeset('RelTol',1e-6,'NonNegative',1:length(y0));
[t,y] = ode45(@typhoid_ode_nepal,1:tspan,y0,options);

%Discard burn-in period
t(1:t0)=[]; 
y(1:t0,:)=[];

D=rep*delta*y(:,al+1:2*al);  
D=sum(D,2);

obs=data;

%%% CALCULATE MODEL FIT (Least-squares or Negative log-likelihood) %%%
%xerr = (obs - D) ;
%LS = xerr'*xerr; 

llikl=zeros(size(obs));
for i=1:size(obs,1)
    for j=1:size(obs,2)
        llikl(i,j)=obs(i,j).*log(D(i,j)) - D(i,j) - sum(log(1:obs(i,j))); %Log-likelihood of data assuming number of cases at time i, age j is Poisson-distributed
    end
end
LL=-sum(sum(llikl)); %Negative log-likelihood

end
