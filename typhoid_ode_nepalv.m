function dydt=typhoid_ode(t,y)
%Differential equations for typhoid model

global B u mu omega r rC alpha delta theta xi gamma epsilon betaw betap q lag al v1 massvacc veff wv tcarr san immig; 

lambdap=betap*(y(al+1:2*al)+r*y(4*al+1:5*al)+rC*y(5*al+1:6*al))/sum(y(1:8*al)); %
if length(betaw)==52
lambdaw=betaw(rem(round(t),52)+1,1)*y(8*al+1)*san(round(t));  %/sum(y(1:6*al));  %
else
lambdaw=betaw(round(t),1)*y(8*al+1)*san(round(t));  %/sum(y(1:6*al));  %
end
for i=1:al
    dydt(i,1)= B(1,i)*(1-v1(round(t),i)) + immig(round(t),i) + omega*epsilon*y(2*al+i) - (lambdap(i)+lambdaw)*y(i) - (u(i)+mu)*y(i); %dS1/dt IF USING BIRTH RATE: log(1+B(round(t),i))/52*sum(St(:))
    dydt(i+al,1)= (lambdap(i)+lambdaw)*(y(i)+(1-veff)*y(i+6*al)) - delta*y(al+i) - (u(i)+mu)*y(al+i); %dI1/dt
    dydt(i+2*al,1)= v1(round(t),i)*B(1,i) + delta*(1-theta)*y(al+i) + delta*(1-theta)*y(4*al+i) - omega*y(2*al+i) + tcarr(round(t))*y(5*al+i) - (u(i)+mu)*y(2*al+i); %dR/dt
    dydt(i+3*al,1)= omega*(1-epsilon)*y(2*al+i) - (lambdap(i)+lambdaw)*y(3*al+i) - (u(i)+mu)*y(3*al+i); %dS2/dt 
    dydt(i+4*al,1)= (lambdap(i)+lambdaw)*(y(3*al+i)+(1-veff)*y(i+7*al)) - delta*y(4*al+i) - (u(i)+mu)*y(4*al+i); %dI2/dt 
    dydt(i+5*al,1)= delta*theta*(y(al+i) + y(4*al+i)) - (tcarr(round(t))+u(i)+mu)*y(5*al+i); %dC/dt
    dydt(i+6*al,1)= -(wv+u(i)+mu)*y(i+6*al); % - omega*epsilon*y(i+1+6*al);
    dydt(i+7*al,1)= -(wv+u(i)+mu)*y(i+7*al); %
    dydt(8*al+1,1)= gamma*sum(y(al+1:2*al)+r*y(4*al+1:5*al)+rC*y(5*al+1:6*al)) - xi*y(8*al+1); %dW/dt
    if i>1 %aging
        dydt(i,1)=dydt(i,1) + u(i-1)*y(i-1)*(1-v1(round(t),i));
        dydt(i+al,1)=dydt(i+al,1) + u(i-1)*y(i+al-1);
        dydt(i+2*al,1)=dydt(i+2*al,1) + u(i-1)*y(i+2*al-1) + v1(round(t),i)*u(i-1)*(y(i-1)+y(i+3*al-1));
        dydt(i+3*al,1)=dydt(i+3*al,1) + u(i-1)*y(i+3*al-1)*(1-v1(round(t),i));
        dydt(i+4*al,1)=dydt(i+4*al,1) + u(i-1)*y(i+4*al-1); 
        dydt(i+5*al,1)=dydt(i+5*al,1) + u(i-1)*y(i+5*al-1);
        %cumulative proportion vaccinated
        dydt(i+6*al,1)= dydt(i+6*al,1) + u(i-1)*y(i+6*al-1); %+y(i+al-1)+y(i+2*al-1)+y(i+3*al-1)+y(i+4*al-1)+y(i+5*al-1)
        dydt(i+7*al,1)= dydt(i+7*al,1) + u(i-1)*y(i+7*al-1);
    end

if massvacc(round(t),i)>0
    dydt(i,1) = dydt(i,1) - massvacc(round(t),i)*y(i);
    dydt(3*al+i,1) = dydt(3*al+i,1) - massvacc(round(t),i)*y(3*al+i);
    dydt(6*al+i,1) = dydt(6*al+i,1) + massvacc(round(t),i)*y(i);
    dydt(7*al+i,1) = dydt(7*al+i,1) + massvacc(round(t),i)*y(3*al+i);
    t=t+1;
end
end