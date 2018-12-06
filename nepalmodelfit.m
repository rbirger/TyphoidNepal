pguess1=[par1.betap(1); par1.betaw; par1.q; par1.epsilon; par1.alpha; par1.rC];
pguess2=[par2.betap(1); par2.betaw; par2.q; par2.epsilon; par2.alpha; par2.rC];
%%
par_nepal1a=fminsearch(@(p) typhoidfit_nepal(p,C_nepal(1:length(rainfall_nepalwk),1),93859,0.0309,rainfall_nepalwk),pguess1)
par_nepal2a=fminsearch(@(p) typhoidfit_nepal(p,C_nepal(1:length(rainfall_nepalwk),1),93859,0.0309,rainfall_nepalwk),pguess2)

%%
[LL1_nepal,Cfit1a_nepal,R01a_nepal]=typhoidfit_nepal(par_nepal1a,C_nepal(1:length(rainfall_nepalwk),1),93859,0.0309,rainfall_nepalwk);
[LL2_nepal,Cfit2a_nepal,R02a_nepal]=typhoidfit_nepal(par_nepal2a,C_nepal(1:length(rainfall_nepalwk),1),93859,0.0309,rainfall_nepalwk);

%%
par_nepal1a=abs(par_nepal1a);
par_nepal2a=abs(par_nepal2a);

pfit1.betap = par_nepal1a(1);
pfit1.betaw = par_nepal1a(2);
pfit1.q = par_nepal1a(3);
pfit1.epsilon = par_nepal1a(4);
pfit1.alpha = par_nepal1a(5);
pfit1.r = par_nepal1a(6);
pfit1.rC = par_nepal1a(6);

pfit2.betap = par_nepal2a(1);
pfit2.betaw = par_nepal2a(2);
pfit2.q = par_nepal2a(3);
pfit2.epsilon = par_nepal2a(4);
pfit2.alpha = par_nepal2a(5);
pfit2.r = par_nepal2a(6);
pfit2.rC = par_nepal2a(6);

%%
figure
hold on
plot(C_nepal(1:length(rainfall_nepalwk),1))
plot(Cfit1a_nepal,'r')
plot(Cfit2a_nepal,'g')

