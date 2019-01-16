function [LS,estim]=seasfit(par,data,vargin)

b0=par(1); a=par(2); phi=par(3);

if isempty(vargin)
    T=12; %tmax; %
else
    T=vargin;
end

tmax=length(data);
estim=zeros(tmax,1);
%data=data';
for t=1:tmax
    %estim(t,1)=b0*(1+a*cos(2*pi*(t-phi*T)/T));
    estim(t,1)=b0+a*cos(2*pi*(t-phi*T)/T);
end

ls=(estim-data).^2;
LS=sum(ls);