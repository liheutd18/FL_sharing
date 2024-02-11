clc 
clear
close all
usernum=10;

%pgsyc=[0.3658 0.3658 0.3658 0.3658 0.3658 0.3658 0.3658 0.8595 0.8595 0.8595 1.3782 1.3782 1.3782 1.3782 1.3782 0.8595 0.8595 0.8595 1.3782 1.3782 1.3782 0.8595 0.8595 0.3658];
%pgbyc=[0.35 0.35 0.35 0.35 0.35 0.35 0.35	0.35 0.35 0.35 0.35	0.35 0.35 0.35 0.35	0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35];
pgsyc=[0.2*ones(1,7) 0.55*ones(1,16)  0.2*ones(1,1) ];
pgbyc=0.1*ones(1,24);
% pgs=1.0*ones(1,24);
% pgb=0.4*ones(1,24);
% perror = -0.05 + 0.1*rand(1,24);
load perror;
perror=0;
pgs=pgsyc+perror;
pgb=pgbyc+perror;

a=zeros(24,24);
SOCmin=0.05;
SOCmax=0.95;
SOCini=0.05;
SOCsd=zeros(25,1);
SDcap=40;
pvbt=0;
SOCsd(1)=SOCini*SDcap;
c=ones(24,1);
k=zeros(24,usernum);
role=zeros(24,usernum);
newload=zeros(24,usernum);
optload=zeros(24,usernum);
cu=zeros(24,usernum);
tu=zeros(24,usernum);
newcu=zeros(24,usernum);
newtu=zeros(24,usernum);
userpro=zeros(1,usernum);
newuserpro=zeros(1,usernum);
bnum=zeros(24,1);
snum=zeros(24,1);
btk=zeros(24,1);
stk=zeros(24,1);
btpv=zeros(24,1);
stpv=zeros(24,1);
price=zeros(24,2);
pro=zeros(24,1);

% [pv,l_0,loadmin,loadmax]=ini();
% pv=1.0*pv/1000;l_0=1.0*l_0/1000;loadmin=loadmin/1000;loadmax=loadmax/1000;
% ntl=l_0-pv;
[pv,l_0,loadmin,loadmax]=Pecan();
ntl=l_0-pv;
% pverror = -2000 + 4000*rand(24,3);
% lerror = -2000 + 4000*rand(24,3);
% nlerror=-1000+2000*rand(24,1);
% load pverror;
% load lerror;

% wcxs=0.9+0.2*rand(24,1);
load wcxs;
% pvyc=pv+pverror;
% lyc=l+lerror;
tl=sum(l_0,2);
tpv=sum(pv,2);
% tlyc=sum(lyc,2);
% tpvyc=sum(pvyc,2);
iniNLt=tl-tpv;


tnl=zeros(24,1);
for i=1:24
tnl(i)=iniNLt(i)*wcxs(i);
end

for i=1:24
    for j=1:i
        a(i,j)=1;
    end
end
b=-a;
% load zj;
A=[a;b];

x0=zeros(24,1); 
b=[SDcap*(SOCmax-SOCini)*c;SDcap*(SOCini-SOCmin)*c]; Aeq=[]; beq=[];
c_rate=2;%1%2
    b=[SDcap*(SOCmax-SOCini)*c;SDcap*(SOCini-SOCmin)*c];
    lb=-(SDcap/c_rate)*c;
    ub=(SDcap/c_rate)*c;

% [Esd,cost]=fmincon(@(x)myfun1(x,tl,tpv),x0,A,b,Aeq,beq,lb,ub);

wb=ones(1,24);
ws=ones(1,24);
oribtpv=zeros(24,1);
oristpv=zeros(24,1);
for T=1:24
    wb(T)=pgs(T)/max(pgs);
    ws(T)=min(pgb)/pgb(T);
        for i=1:usernum
            if(l_0(T,i)>=pv(T,i))
                k(T,i)=(l_0(T,i)+1)*pgs(T);
                role(T,i)=1;
                bnum(T)=bnum(T)+1;
                oribtpv(T)=oribtpv(T)+l_0(T,i)-pv(T,i);
                btk(T)=btk(T)+k(T,i);
            else
                k(T,i)=(l_0(T,i)+1)*pgb(T);
                role(T,i)=0;
                snum(T)=snum(T)+1;
                oristpv(T)=oristpv(T)+l_0(T,i)-pv(T,i);
                stk(T)=stk(T)+k(T,i);
            end
        end
end
oristpv=-oristpv;
wb=[0.2*ones(1,7) 0.55*ones(1,16)  0.2*ones(1,1) ];
ws=0.1*ones(1,24);
%tnl=tnl.*(0.8+0.2*rand(24,1));
tic;
options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',10000);
nonlcon=[];
[Esd,cost]=fmincon(@(x)myfuntest(x,tnl,wb,ws),x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
toc;
eff=1;
cost1=wb(1)*max(tnl(1)+Esd(1)/eff,0)+ws(1)*min(tnl(1)+Esd(1)*eff,0);
DDESG=Esd+tnl;

for i=2:25
    SOCsd(i)=SOCsd(i-1)+Esd(i-1);
end

    gamma = 0.5;    % discount factor  % TODO : we need learning rate schedule
    alpha = 0.5;    % learning rate    % TODO : we need exploration rate schedule
    epsilon = 0.2;  % exploration probability (1-epsilon = exploit / epsilon = explore)    
tic;
    J=5000;
    Q = zeros(1,J);
    Q1= -200*ones(24,J);
    Q_all=0;
    l=l_0;
for j=1:J      
    lihe=1;
    %l=l_0;
    disp(j)
    for T=1:24
%        clearvars Q1

   % x0=[pgs(T);pgb(T)]; A=[]; b=[]; Aeq=[]; beq=[]; 
    lb=pgb(T)*ones(2,1); ub=pgs(T)*ones(2,1);  
    pb_action=[lb(1):0.01:ub(1)];
    ps_action=[lb(2):0.01:ub(2)];
    action=1:1:length(pb_action)*length(ps_action); 
    ActionMAT= cell(length(pb_action),length(ps_action));
    for m=1:length(pb_action)
        for n=1:length(ps_action)
            ActionMAT(m,n)={[pb_action(1,m),ps_action(1,n)]};
        end
    end    
%    epsilon=1-1/(J+1-j);
    %epsilon=1/J;
%    for o=1:length(action)
        r=rand; % get 1 uniform random number
        x=sum(r>=cumsum([0, 1-epsilon, epsilon])); % check it to be in which probability area
        %choose either explore or exploit
        if x == 1   % exploit;
                [C,I]=max(Q1(T,:),[],2);    
            if I<j                
                bestaction=Pricerecordall{T,I};           
                row=bestaction(1)*100-9;col=bestaction(2)*100-9;
                action_idx=(row-1)*length(ActionMAT)+col; %max;
            else
                action_idx=datasample([1:1:length(action)],1);
            end
        else        % explore
            action_idx=datasample([1:1:length(action)],1);
%            [r_mat,c_mat]=find(IndexMAT==I);  
%             Paction=ActionMAT{r_mat,c_mat}; 
%             choose 1 action randomly (uniform random distribution)
%            action_idx=I;
        end    
        
%    action_idx = find(action==o); % id of the chosen action
    % observe the next state and next reward ** there is no reward matrix
    
    IndexMAT=[1:1:length(pb_action)*length(ps_action)];    IndexMAT=reshape(IndexMAT,length(ps_action),length(pb_action))';
    [r_mat,c_mat]=find(abs(IndexMAT-action_idx)<0.01);% [i j]=find(abs(a-0.6678)<1e-5) 
    Paction=ActionMAT{r_mat,c_mat};    
    price=Paction;
        if bnum(T)==usernum
            price(1)=ub(1);price(2)=lb(2);
%         elseif snum(T)==usernum
%             price(1)=lb(1);price(2)=lb(2);
        elseif price(1)<price(2)
            price(1)=price(2);
        end
% % demand response 
        for i=1:usernum
            if(role(T,i)==1)
                optload(T,i)=k(T,i)/price(1)-1;
            else
                optload(T,i)=k(T,i)/price(2)-1;
            end
            if(optload(T,i)>loadmax(T,i))
                newload(T,i)=loadmax(T,i);
            elseif(optload(T,i)<loadmin(T,i))
                newload(T,i)=loadmin(T,i);
            else
                newload(T,i)=optload(T,i);
            end
        end           

        [Q1(T,j),r]=myfun2future(newload(T,:),l(T:end,:),pv(T:end,:),price,stk(T:end),btk(T:end),snum(T:end),bnum(T:end),stpv(T:end),btpv(T:end),pgb(T:end),pgs(T:end),Esd(T:end));   
    
%    end
%             [C,I]=max(Q1(T,:),[],2); 
% 
%             [r_mat,c_mat]=find(IndexMAT==I);
%     %        reward=r(I);
%             Paction=ActionMAT{r_mat,c_mat};    
%             price=Paction;

%         r=rand; % get 1 uniform random number
%         x=sum(r>=cumsum([0, 1-epsilon, epsilon])); % check it to be in which probability area
%         %choose either explore or exploit
%         if x == 1   % exploit
%             price=Paction;
%         else        % explore
%            I=datasample([1:1:length(Q1)],1);
%            [r_mat,c_mat]=find(IndexMAT==I);  
%             Paction=ActionMAT{r_mat,c_mat}; 
%             % choose 1 action randomly (uniform random distribution)
%         end
        for i=1:usernum
            if(role(T,i)==1)
                optload(T,i)=k(T,i)/price(1)-1;
            else
                optload(T,i)=k(T,i)/price(2)-1;
            end
            if(optload(T,i)>loadmax(T,i))
                newload(T,i)=loadmax(T,i);
            elseif(optload(T,i)<loadmin(T,i))
                newload(T,i)=loadmin(T,i);
            else
                newload(T,i)=optload(T,i);
            end
        end
        % reward
        bamount=0;samount=0;        
        for i=1:usernum
            if(newload(T,i)>0)
                role(T,i)=1;
                bamount=bamount+newload(T,i);
            else
                role(T,i)=0;
                samount=samount+newload(T,i);
            end
        end
        
        dr=l(T,:)-newload(T,:);
        if T<24
            avedr=dr/(24-T);
            l(T+1:end,:)=l(T+1:end,:)-avedr;            
        else
            avedr=0;
        end

        l(T,:)=newload(T,:);
        Predntl=sum(l,2)-sum(pv,2);
%          load lstm0216.mat %net sig mu
%          net = resetState(net);
%          numTimeStepsTest = 24-T;
%          for i = 1:numTimeStepsTest
%          [net,YPred(:,i)] = predictAndUpdateState(net,Predntl(:,i),'ExecutionEnvironment','cpu');
%          end
% UPDATE WITH LSTM
%         XTrain=TotalNTL(1:24*d)';
%         XTrain=Predntl;
%         net = predictAndUpdateState(net,XTrain);
%         YPred = [];
%         numTimeStepsTest = 24-T;
%         for i = 1:numTimeStepsTest
%             [net,YPred(:,i)] = predictAndUpdateState(net,Predntl(:,i),'ExecutionEnvironment','cpu');
%         end        
%        YPred = sig*YPred + mu;        
        
        DeltaE=bamount+samount+Esd(T);
        %reward=r;
        reward=price(1)*bamount+price(2)*samount-pgs(T)*max(DeltaE,0)-pgb(T)*min(DeltaE,0)-0.15*(max(Esd(T),0));
        %[Q1_opt,quantity_opt,l_opt]=myfun2future(l(T:end,:),pv(T:end,:),price,stk(T:end),btk(T:end),snum(T:end),bnum(T:end),stpv(T:end),btpv(T:end),pgb(T:end),pgs(T:end),Esd(T:end));       
        %reward=stk(T)+btk(T)-price(1).*(bnum(T)+btpv(T))-price(2).*(snum(T)+stpv(T))...
%            +pgb(T)'.*max((stpv(T)+btpv(T)+snum(T)+bnum(T)-(btk(T)./price(1)))-(stk(T)./price(2)-Esd(T)),0)+pgs(T)'.*min((stpv(T)+btpv(T)+snum(T)+bnum(T)-(btk(T)./price(1))-(stk(T)./price(2))-Esd(T)),0);%-0.15*sum(max(Esd,0));
        %reward=price(1)*quantity_opt(1,1)+price(2)*quantity_opt(1,2)-pgs(T)*max(quantity_opt(1,3),0)-pgb(T)*min(quantity_opt(1,3),0);
        
        Pricerecordall{T,j}=price;
%         reward=max(Q1);   
%         Q(j)=Q(j)+theta*(reward+gamma*max(Q1)-Q(j));
         
        %load update
%         if T<23
%             l(T+1:end,:)=l_opt(2:end,:);
%         elseif T==24
%             l(24,:)=l_opt(end,:);
%         end      
%         for i=1:usernum
% %             if(role(T,i)==1)
% %                 l(T,i)=k(T,i)/price(1)-1;
% %             else
% %                 l(T,i)=k(T,i)/price(2)-1;
% %             end
% 
%             if(l(T,i)>loadmax(T,i))
%                 l(T,i)=loadmax(T,i);
%             elseif(l(T,i)<loadmin(T,i))
%                 l(T,i)=loadmin(T,i);
%             else 
%                 l(T,i)=l(T,i);
%             end
%         end    
    
    theta=0.1;gamma=0.95 ;
    Q(j)=Q(j)+theta*(reward+gamma*max(Q1(T,:))-Q(j));
    Q_all=[Q_all;Q(j)];
    end
    %l=newload;
    li=1;
end

plot(Q);
toc;
price_final=cell2mat(Pricerecordall(:,4997));
newNLt=sum(newload,2)-tpv+Esd;
%newNLt(14)=-35;% why -35
%optNLt=sum(optload,2)-tpv+Esd;

optEbuy=zeros(24,1);
optEsell=zeros(24,1);
for T=1:24
    for i=1:usernum
        if(newload(T,i)-pv(T,i)>=0)
            optEsell(T)=optEsell(T)+newload(T,i)-pv(T,i);
        else
            optEbuy(T)=optEbuy(T)-(newload(T,i)-pv(T,i));
        end
    end
end

for T=1:24
    if(newNLt(T)<0)
        optMGOpro(T)=price_final(T,1)*optEsell(T)-price_final(T,2)*optEbuy(T)+pgb(T)*(-newNLt(T))-0.15*(max(Esd(T),0));
    else
        optMGOpro(T)=price_final(T,1)*optEsell(T)-price_final(T,2)*optEbuy(T)+pgs(T)*(-newNLt(T))-0.15*(max(Esd(T),0));
    end
end


for T=1:24
    for i=1:usernum
        userpro(T,i)=k(T,i)*log(1+l_0(T,i))-pgs(T)*max(l_0(T,i)-pv(T,i),0)-pgb(T)*min(l_0(T,i)-pv(T,i),0)+pvbt*pv(T,i);
        cu(T,i)=k(T,i)*log(1+l(T,i));
        tu(T,i)=-pgs(T)*max(l(T,i)-pv(T,i),0)-pgb(T)*min(l(T,i)-pv(T,i),0)+pvbt*pv(T,i);
        newuserpro(T,i)=k(T,i)*log(1+newload(T,i))-price_final(T,1)*max(newload(T,i)-pv(T,i),0)-price_final(T,2)*min(newload(T,i)-pv(T,i),0)+pvbt*pv(T,i);
        newcu(T,i)=k(T,i)*log(1+newload(T,i));
        newtu(T,i)=-price_final(T,1)*max(newload(T,i)-pv(T,i),0)-price_final(T,2)*min(newload(T,i)-pv(T,i),0)+pvbt*pv(T,i);
    end
end
% Table results
Final_utilitynew=sum(newuserpro,1);Final_utilityold=sum(userpro,1);
Final2uti=[Final_utilityold;Final_utilitynew];
Final_profit=sum(optMGOpro);
for i=1:usernum
    NL(:,i)=l_0(:,i)-pv(:,i);
end
Final_L0=sum(NL,2);
for i=1:usernum
    newNL(:,i)=newload(:,i)-pv(:,i);
end
Final_Lnew=sum(newNL,2);
Final_LnewwithES=sum(Final_Lnew+Esd,2);
figure(1)
plot(Final_L0)
hold on
plot(Final_Lnew)
hold on
plot(Final_LnewwithES)
% realEbuy=zeros(24,1);
% realEsell=zeros(24,1);
% for T=1:24
%     for i=1:usernum
%         if(newload(T,i)-pv(T,i)>=0)
%             realEsell(T)=realEsell(T)+newload(T,i)-pv(T,i);
%         else
%             realEbuy(T)=realEbuy(T)-(newload(T,i)-pv(T,i));
%         end
%     end
% end

% realMGOpro=zeros(24,1);
% for T=1:24
%     if(realEsell(T)>realEbuy(T))
%         for i=1:usernum
%             if(newload(T,i)-pv(T,i)<0)
%                 if(newload(T,i)<optload(T,i))
%                     realuserpro(T,i)=k(T,i)*log(1+newload(T,i))+price(T,2)*(pv(T,i)-newload(T,i))+pvbt*pv(T,i);
%                     realMGOpro(T)=realMGOpro(T)-price(T,2)*(pv(T,i)-newload(T,i));
%                 else
%                     realuserpro(T,i)=k(T,i)*log(1+newload(T,i))+price(T,2)*(pv(T,i)-newload(T,i))+(price(T,2)-pgb(T))*(optload(T,i)-newload(T,i))+pvbt*pv(T,i);
%                     realMGOpro(T)=realMGOpro(T)-(price(T,2)*(pv(T,i)-newload(T,i))+(price(T,2)-pgb(T))*(optload(T,i)-newload(T,i)));
%                 end
%             else
%                 if(newload(T,i)>optload(T,i))
%                     realuserpro(T,i)=k(T,i)*log(1+newload(T,i))-price(T,1)*(optload(T,i)-pv(T,i))-pgs(T)*(newload(T,i)-optload(T,i))+pvbt*pv(T,i);
%                     realMGOpro(T)=realMGOpro(T)+(price(T,1)*(optload(T,i)-pv(T,i))+pgs(T)*(newload(T,i)-optload(T,i)));
%                 else
%                     realuserpro(T,i)=k(T,i)*log(1+newload(T,i))-price(T,1)*(newload(T,i)-pv(T,i))+pvbt*pv(T,i);
%                     realMGOpro(T)=realMGOpro(T)+price(T,1)*(newload(T,i)-pv(T,i));
%                 end
%             end
%         end
%     else
%         for i=1:usernum
%             if(newload(T,i)-pv(T,i)<0)
%                 if(newload(T,i)<optload(T,i))
%                     realuserpro(T,i)=k(T,i)*log(1+newload(T,i))+price(T,2)*(pv(T,i)-optload(T,i))+pmb(T)*(optload(T,i)-newload(T,i))+pvbt*pv(T,i);
%                     realMGOpro(T)=realMGOpro(T)-(price(T,2)*(pv(T,i)-optload(T,i))+pmb(T)*(optload(T,i)-newload(T,i)));
%                 else
%                     realuserpro(T,i)=k(T,i)*log(1+newload(T,i))+price(T,2)*(pv(T,i)-newload(T,i))+pvbt*pv(T,i);
%                     realMGOpro(T)=realMGOpro(T)-price(T,2)*(pv(T,i)-newload(T,i));
%                 end
%             else
%                 if(newload(T,i)>optload(T,i))
%                     realuserpro(T,i)=k(T,i)*log(1+newload(T,i))-price(T,1)*(newload(T,i)-pv(T,i))+pvbt*pv(T,i);
%                     realMGOpro(T)=realMGOpro(T)+price(T,1)*(newload(T,i)-pv(T,i));
%                 else
%                     realuserpro(T,i)=k(T,i)*log(1+newload(T,i))-price(T,1)*(newload(T,i)-pv(T,i))+(pgs(T)-price(T,1))*(newload(T,i)-optload(T,i))+pvbt*pv(T,i);
%                     realMGOpro(T)=realMGOpro(T)+(price(T,1)*(newload(T,i)-pv(T,i))-(pgs(T)-price(T,1))*(newload(T,i)-optload(T,i)));
%                 end
%             end
%         end  
%     end
%     if(-newNLt(T)-Esd(T)>0)
%         realMGOpro(T)=realMGOpro(T)+pgb(T)*(-newNLt(T));
%     else
%         realMGOpro(T)=realMGOpro(T)+pgs(T)*(-newNLt(T));
%     end
% end       

% figure(1)
% plot(1:1:24,newuserpro-userpro);
% % figure(2)
% % plot(1:1:24,load,1:1:24,newload);
% % figure(3)
% plot(1:1:24,iniNLt,1:1:24,newNLt);

% figure(2)
% plot(1:1:24,newNL-NL);