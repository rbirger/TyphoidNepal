function dydt=typhoid_ode_nepal(t,y)
%Differential equations for typhoid model

global B u mu omega r alpha delta theta xi gamma epsilon betaw betap q lag al v1 massvacc immig san t1 t2 rt; % noise; rC 

rC=r;
theta2=theta; %zeros(1,al); %

trans=1;
%dur=4;
if round(t)>=t1 && round(t)<t2
    trans=1+(rt-1)*(round(t)-t1)/(t2-t1);
    %delta=1/(dur*trans);
elseif round(t)>=t2
    trans=rt;
    %delta=1/(dur*rt);
end

lambdap=betap*(y(al+1:2*al)+r*y(4*al+1:5*al)+rC*y(5*al+1:6*al))/sum(y(1:6*al))*trans; %(1+q*cos(2*pi*(t-lag)/52.18))* noise(round(t))*
lambdaw=max(0,(1+q*cos(2*pi*(t-lag)/52.18))*betaw(round(t),1)*(y(6*al+1))*san(round(t),1)*trans)/sum(y(1:6*al)); %noise(round(t))*
%if length(betaw)==52
%lambdaw=betaw(rem(round(t),52)+1,1)*y(6*al+1);  %/sum(y(1:6*al));  %
%else
%lambdaw=betaw(round(t),1)*y(6*al+1)*san(round(t),1)*trans; %;  %/sum(y(1:6*al));  %
%end

for i=1:al
    dydt(i,1)= B(1,i)*(1-v1(round(t),i)) + immig(round(t),i) + omega*epsilon*y(2*al+i) - (lambdap(i)+lambdaw)*y(i) - (u(i)+mu)*y(i); %dS1/dt IF USING BIRTH RATE: log(1+B(round(t),i))/52*sum(St(:))
    dydt(i+al,1)= (lambdap(i)+lambdaw)*y(i) - delta*y(al+i) - (u(i)+mu)*y(al+i); %dI1/dt
    dydt(i+2*al,1)= B(1,i)*v1(round(t),i) + delta*(1-alpha-theta(i))*y(al+i) + delta*(1-theta2(i))*y(4*al+i) - omega*y(2*al+i) - (u(i)+mu)*y(2*al+i); %dR/dt
    dydt(i+3*al,1)= omega*(1-epsilon)*y(2*al+i) - (lambdap(i)+lambdaw)*y(3*al+i) - (u(i)+mu)*y(3*al+i); %dS2/dt 
    dydt(i+4*al,1)= (lambdap(i)+lambdaw)*y(3*al+i) - delta*y(4*al+i) - (u(i)+mu)*y(4*al+i); %dI2/dt 
    dydt(i+5*al,1)= delta*(theta(i)*y(al+i) + theta2(i)*y(4*al+i)) - (u(i)+mu)*y(5*al+i); %dC/dt
    dydt(6*al+1,1)= gamma*sum(y(al+1:2*al)+r*y(4*al+1:5*al)+rC*y(5*al+1:6*al)) - xi*y(6*al+1); %dW/dt
    dydt(i+6*al+1,1)= v1(round(t),i)*B(1,i) - u(i)*y(i+1+6*al); % - omega*epsilon*y(i+1+6*al);
    if i>1 %aging
        dydt(i,1)=dydt(i,1) + u(i-1)*y(i-1)*(1-v1(round(t),i));
        dydt(i+al,1)=dydt(i+al,1) + u(i-1)*y(i+al-1);
        dydt(i+2*al,1)=dydt(i+2*al,1) + u(i-1)*y(i+2*al-1) + u(i-1)*(y(i-1) + y(i+3*al-1))*v1(round(t),i);
        dydt(i+3*al,1)=dydt(i+3*al,1) + u(i-1)*y(i+3*al-1)*(1-v1(round(t),i));
        dydt(i+4*al,1)=dydt(i+4*al,1) + u(i-1)*y(i+4*al-1); 
        dydt(i+5*al,1)=dydt(i+5*al,1) + u(i-1)*y(i+5*al-1);
        %cumulative proportion vaccinated
        dydt(i+1+6*al,1)= dydt(i+1+6*al,1) + v1(round(t),i)*u(i-1)*(y(i-1) + y(i+3*al-1)) + u(i-1)*y(i+1+6*al-1); %+y(i+al-1)+y(i+2*al-1)+y(i+3*al-1)+y(i+4*al-1)+y(i+5*al-1)
    end
end
if massvacc(round(t),1)>0
    dydt(1:al,1) = dydt(1:al,1) - massvacc(round(t),1)*y(1:al);
    dydt(2*al+1:3*al,1) = dydt(2*al+1:3*al,1) + massvacc(round(t),1)*(y(1:al) + y(3*al+1:4*al));
    dydt(3*al+1:4*al,1) = dydt(3*al+1:4*al,1) - massvacc(round(t),1)*y(3*al+1:4*al);
    dydt(6*al+2:7*al+1,1) = dydt(6*al+2:7*al+1,1) + massvacc(round(t),1)*(y(1:al) + y(3*al+1:4*al));
    t=t+1;
end