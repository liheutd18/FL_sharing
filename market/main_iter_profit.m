clc 
clear
close all
load('resultsday313.mat')
clearvars -except estimatepv

usernum=10;
tic

pgsyc=[0.3658 0.3658 0.3658 0.3658 0.3658 0.3658 0.3658 0.8595 0.8595 0.8595 1.3782 1.3782 1.3782 1.3782 1.3782 0.8595 0.8595 0.8595 1.3782 1.3782 1.3782 0.8595 0.8595 0.3658];
pgbyc=[0.35 0.35 0.35 0.35 0.35 0.35 0.35	0.35 0.35 0.35 0.35	0.35 0.35 0.35 0.35	0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35];
% pgs=1.0*ones(1,24);
% pgb=0.4*ones(1,24);
% perror = -0.05 + 0.1*rand(1,24);
load perror;
perror=zeros(1,24);
pgs=pgsyc+perror;
pgb=pgbyc+perror;

a=zeros(24,24);
SOCmin=0.05;
SOCmax=0.95;
SOCini=0.05;
SOCsd=zeros(25,1);
SDcap=40*1000;
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
bk=zeros(24,usernum);
sk=zeros(24,usernum);
btv=zeros(24,usernum);
stv=zeros(24,usernum);
btpv=zeros(24,1);
stpv=zeros(24,1);
price=zeros(24,2);
priceb=zeros(24,usernum);
prices=zeros(24,usernum);
finalprices=zeros(24,usernum);
pro=zeros(24,1);
eff=0.95;

     [ntl l pv] = Pecan();
     ntl = ntl(242:242+23,:);
     l = l(242:242+23,:);
     pv = pv(242:242+23,:); 
     pv_real = pv;
     ntl_real = ntl;
     
     Data1 =readtable('Forecastingresults2lists.xlsx','Sheet','Agg');
     Data2=readtable('Forecastingresults2lists.xlsx','Sheet','Individual');     
     AGGNetLoad=table2array(Data1(:,2:25)); 
     AGGNetLoad = AGGNetLoad'; AGGNetLoad = AGGNetLoad*1000;
     INDNetLoad=table2array(Data2(:,3:26)); 
     INDNetLoad = INDNetLoad'; INDNetLoad = INDNetLoad*1000;

     
     pv_dis = estimatepv(242:242+23,:);
     pv_dis(pv_dis<0)=0; pv_dis=pv_dis*1000;
     pv_real(pv_real<0)=0;
     pv_real=pv_real*1000;      pv=pv*1000; 
     ntl_real=ntl_real*1000; 
     l_real = ntl_real + pv_real;   
     
%%     
%     for forecasting_selection1 = 1:3 % 2 or 3
     forecasting_selection1 = 1;
     ntl_dis = [];
     for foreseq = 1 : 10
        ntl_dis = [ntl_dis INDNetLoad(:,(foreseq-1)*3 + forecasting_selection1)];
     end
     
     for forecasting_selection2 = 1:5 % 2 or 3
     
%     forecasting_selection2 = 1; % 2 or 3 4 5
     tnl = AGGNetLoad (:,forecasting_selection2);
     sum(ntl_dis,2) - tnl

     
%%%% BTM MAPE     
     error = mean(abs(pv_real-pv_dis)./(pv_real+1)); error(5) = 0; % user 5 has no PV
     
     l=ntl_dis + pv_dis;

     loadmin=0.8*l_real;
     loadmax=1.2*l_real;


% tnl=AGGNetLoad_predDA';tnl=tnl*1000;
% tnl=sum(ntl,2);

     tl=sum(l_real,2);
     tpv=sum(pv_real,2);

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


wb=ones(1,24);
ws=ones(1,24);
for T=1:24
    wb(T)=pgs(T);%/max(pgs);
    ws(T)=min(pgb);%/pgb(T);
end
options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',10000);
nonlcon=[];
tic;

        c_rate=2; %1%2
        b=[SDcap*(SOCmax-SOCini)*c;SDcap*(SOCini-SOCmin)*c];
        lb=-(SDcap/c_rate)*c;
        ub=(SDcap/c_rate)*c;
        [Esd,cost,cost1]=fmincon(@(x)myfuntest(x,tnl,wb,ws),x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
 
toc;
New_NL=max(tnl+Esd/eff,0)+min(tnl+Esd*eff,0);
cost_all5(forecasting_selection2) = cost;

for i=2:25
    SOCsd(i)=SOCsd(i-1)+Esd(i-1); 
end

% f=wb*max(Esd,0)+ws*min(Esd,0);
% Esd(Esd<0)=0;
% Esd_fee=sum(sum(Esd.*(wb-ws)))/10000; 
% cost(1:7)=Esd(1:7)*0.3658/eff
% cost(8:10)=Esd(8:10)*0.8595*eff
% cost(11)=Esd(11)*1.3782*eff
% cost(12:18)=Esd(12:18)*0.35/eff
% cost(19:21)=Esd(19:21)*1.3782*eff
% cost(22:23)=Esd(22:23)*0.8595*eff
% cost(24)=Esd(24)*0.3658*eff
% Cost_all=sum(cost);

for T=1:24
    for i=1:usernum
        if(l(T,i)>=pv_dis(T,i))
            k(T,i)=(l(T,i)+1)*pgs(T);
            role(T,i)=1;
            bnum(T)=bnum(T)+1;
            btpv(T)=btpv(T)+pv_dis(T,i);
            btk(T)=btk(T)+k(T,i);
            bk(T,i)=k(T,i);
            btv(T,i)=pv_dis(T,i);
        else
            k(T,i)=(l(T,i)+1)*pgb(T);
            role(T,i)=0;
            snum(T)=snum(T)+1;
            stpv(T)=stpv(T)+pv_dis(T,i);
            stk(T)=stk(T)+k(T,i);
            sk(T,i)=k(T,i);
            stv(T,i)=pv_dis(T,i);
        end
    end
    %%kreal
    for i=1:usernum
        if(l_real(T,i)>=pv(T,i))
            kreal(T,i)=(l_real(T,i)+1)*pgs(T);
             role_real(T,i)=1;
%             bnum(T)=bnum(T)+1;
%             btpv(T)=btpv(T)+pv(T,i);
%             btk(T)=btk(T)+kreal(T,i);
%             bk(T,i)=kreal(T,i);
%             btv(T,i)=pv(T,i);
        else
            kreal(T,i)=(l_real(T,i)+1)*pgb(T);
             role_real(T,i)=0;
%             snum(T)=snum(T)+1;
%             stpv(T)=stpv(T)+pv(T,i);
%             stk(T)=stk(T)+kreal(T,i);
%             sk(T,i)=kreal(T,i);
%             stv(T,i)=pv(T,i);
        end
    end
    
    kall=[k(T,:);bk(T,:);sk(T,:)];
    bstv=[btv(T,:);stv(T,:)];
    x0=[pgs(T)*ones(1,usernum);pgb(T)*ones(1,usernum)]; 
    A=[]; b=[]; Aeq=[]; beq=[]; 
    lb=pgb(T)*ones(2,usernum); ub=pgs(T)*ones(2,usernum);
    
    [price1,pro1,exitflag]=fmincon(@(x)myfun2_20(x,role(T,:),stk(T),btk(T),kall,bstv,snum(T),bnum(T),stpv(T),btpv(T),pgb(T),Esd(T)),x0,A,b,Aeq,beq,lb,ub,@(x)mycon1_20(x,btk(T),stk(T),kall,snum(T),bnum(T),stpv(T),btpv(T),Esd(T)));
    Flag(T,1)=exitflag;
    [price2,pro2,exitflag]=fmincon(@(x)myfun3_20(x,role(T,:),stk(T),btk(T),kall,bstv,snum(T),bnum(T),stpv(T),btpv(T),pgs(T),Esd(T)),x0,A,b,Aeq,beq,lb,ub,@(x)mycon2_20(x,btk(T),stk(T),kall,snum(T),bnum(T),stpv(T),btpv(T),Esd(T)));
    Flag(T,2)=exitflag;
    
    %sum(kall(2,:)./x(1,:))+sum(kall(3,:)./x(2,:))-snum-bnum-stpv-btpv+Esd
        %if(btk(T)/price1(1)+stk(T)/price1(2)-snum(T)-bnum(T)-stpv(T)-btpv(T)+Esd(T)>=0)
         if sum((kall(2,:)./price1(1,:)+kall(3,:)./price2(1,:))-snum(T)-bnum(T)-stpv(T)-btpv(T)+Esd(T)>=0)   
            pro(T)=pro2;
            priceb(T,:)=price2(1,:);
            prices(T,:)=price2(2,:);
        %elseif(btk(T)/price2(1)+stk(T)/price2(2)-snum(T)-bnum(T)-stpv(T)-btpv(T)+Esd(T)<=0)
         elseif sum((kall(2,:)./price2(1,:)+kall(3,:)./price1(1,:))-snum(T)-bnum(T)-stpv(T)-btpv(T)+Esd(T)>=0)    
            pro(T)=pro1;
            priceb(T,:)=price1(1,:);
            prices(T,:)=price1(2,:);
        else
            if(pro1<pro2)
                pro(T)=pro1;
                priceb(T,:)=price1(1,:);
                prices(T,:)=price1(2,:);
            else
                pro(T)=pro2;
                priceb(T,:)=price2(1,:);
                prices(T,:)=price2(2,:);
            end
         end
        
    %% need real K    
    for i=1:usernum
        if(role_real(T,i)==1)
%            if (role(T,i)==1)
                optload(T,i)=kreal(T,i)/priceb(T,i)-1;
                finalprices(T,i)=priceb(T,i);
%             else
%                 optload(T,i)=kreal(T,i)/ws(T)-1;
%                 finalprices(T,i)=-ws(T);                
%             end            
        else
%            if (role(T,i)==0)
                optload(T,i)=kreal(T,i)/prices(T,i)-1;
                finalprices(T,i)=-prices(T,i);
%             else
%                 optload(T,i)=kreal(T,i)/wb(T)-1;
%                 finalprices(T,i)=wb(T);                
%             end               
        end
        
        if(optload(T,i)>loadmax(T,i))
            newload(T,i)=loadmax(T,i);
        elseif(optload(T,i)<loadmin(T,i))
            newload(T,i)=loadmin(T,i);
        else
            newload(T,i)=optload(T,i);
        end
    end
    
end

newNLt=sum(newload,2)-tpv+Esd;
optNLt=sum(optload,2)-tpv+Esd;

optEbuy=zeros(24,1);
optEsell=zeros(24,1);
for T=1:24
    for i=1:usernum
        optntl(T,i)=optload(T,i)-pv(T,i);
        if(optload(T,i)-pv(T,i)>=0)
            optEsell(T)=optEsell(T)+optload(T,i)-pv(T,i);
        else
            optEbuy(T)=optEbuy(T)-(optload(T,i)-pv(T,i));
        end
    end
end

%%

for T=1:24
    if(+newload(T)+Esd(T)<0)
        %optMGOpro(T)=price(T,1)*optEsell(T)-price(T,2)*optEbuy(T)+pgb(T)*(-optNLt(T));
        optMGOpro(T)=sum((newload(T,:)-pv(T,:)).*finalprices(T,:))+pgb(T)*(-newNLt(T)-Esd(T));
    else
        %optMGOpro(T)=price(T,1)*optEsell(T)-price(T,2)*optEbuy(T)+pgs(T)*(-optNLt(T));
        optMGOpro(T)=sum((newload(T,:)-pv(T,:)).*finalprices(T,:))+pgs(T)*(-newNLt(T)-Esd(T));
    end
end


for T=1:24
    for i=1:usernum
        userpro(T,i)=kreal(T,i)*log(1+l_real(T,i))-pgs(T)*max(l_real(T,i)-pv_real(T,i),0)-pgb(T)*min(l_real(T,i)-pv_real(T,i),0)+pvbt*pv_real(T,i);
        cu(T,i)=kreal(T,i)*log(1+l_real(T,i));
        tu(T,i)=-pgs(T)*max(l_real(T,i)-pv_real(T,i),0)-pgb(T)*min(l_real(T,i)-pv_real(T,i),0)+pvbt*pv_real(T,i);
        
        newuserpro(T,i)=kreal(T,i)*log(1+newload(T,i))-finalprices(T,i)*max(newload(T,i)-pv_real(T,i),0)+finalprices(T,i)*min(newload(T,i)-pv_real(T,i),0)+pvbt*pv_real(T,i);
%         newcu(T,i)=k(T,i)*log(1+newload(T,i));
%         newtu(T,i)=-finalprices(T,i)*max(newload(T,i)-pv(T,i),0)-finalprices(T,2)*min(newload(T,i)-pv(T,i),0)+pvbt*pv(T,i);
    end
end



plot(priceb)
hold on
plot(prices)
plot(finalprices)
plot(1:1:24,newuserpro-userpro);
pro_agent=sum(optMGOpro)/10000;

benefit=sum(newuserpro-userpro,1);%./sum(userpro,1);
benefit_per=100*sum(newuserpro-userpro,1)./sum(userpro,1);
uti1=sum(newuserpro,1)/10000;
uti2=sum(cu+tu,1)/10000;

Benefit_com(forecasting_selection2,:) =  [pro_agent benefit_per];
     end
toc