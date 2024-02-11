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
     
     forecasting_selection1 = 1; % 2 or 3
     ntl_dis = [];
     for foreseq = 1 : 10
        ntl_dis = [ntl_dis INDNetLoad(:,(foreseq-1)*3 + forecasting_selection1)];
     end
     forecasting_selection2 = 1; % 2 or 3 4 5
     tnl = AGGNetLoad (:,forecasting_selection2);
    
     pv_real(pv_real<0)=0;
     pv_real=pv_real*1000;      pv=pv*1000; 
     ntl_real=ntl_real*1000; 
     l_real = ntl_real + pv_real;
     
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
%     wb=[0.2*ones(1,7) 0.55*ones(1,16)  0.2*ones(1,1) ];
%     ws=0.1*ones(1,24);
    c_rate=2;%1%2
%    for i=1:500
%        SDcap=40*1000;
        b=[SDcap*(SOCmax-SOCini)*c;SDcap*(SOCini-SOCmin)*c];
        lb=-(SDcap/c_rate)*c;
        ub=(SDcap/c_rate)*c;
        [Esd,cost,cost1]=fmincon(@(x)myfuntest(x,tnl,wb,ws),x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
        cost=cost1+0.15*i;
%    end
toc;
New_NL=max(tnl+Esd/eff,0)+min(tnl+Esd*eff,0);

for i=2:25
    SOCsd(i)=SOCsd(i-1)+Esd(i-1);
end

for T=1:24
    for i=1:usernum
        if(l(T,i)>=pv(T,i))
            k(T,i)=(l(T,i)+1)*pgs(T);
            role(T,i)=1;
            bnum(T)=bnum(T)+1;
            btpv(T)=btpv(T)+pv(T,i);
            btk(T)=btk(T)+k(T,i);
            bk(T,i)=k(T,i);
            btv(T,i)=pv(T,i);
        else
            k(T,i)=(l(T,i)+1)*pgb(T);
            role(T,i)=0;
            snum(T)=snum(T)+1;
            stpv(T)=stpv(T)+pv(T,i);
            stk(T)=stk(T)+k(T,i);
            sk(T,i)=k(T,i);
            stv(T,i)=pv(T,i);
        end
    end    
    x0=[pgs(T);pgb(T)]; A=[]; b=[]; Aeq=[]; beq=[]; lb=pgb(T)*ones(2,1); ub=pgs(T)*ones(2,1);
    
    [price1,pro1,exitflag]=fmincon(@(x)myfun2(x,stk(T),btk(T),snum(T),bnum(T),stpv(T),btpv(T),pgb(T),Esd(T)),x0,A,b,Aeq,beq,lb,ub,@(x)mycon1(x,btk(T),stk(T),snum(T),bnum(T),stpv(T),btpv(T),Esd(T)));
    Flag(T,1)=exitflag;
    [price2,pro2,exitflag]=fmincon(@(x)myfun3(x,stk(T),btk(T),snum(T),bnum(T),stpv(T),btpv(T),pgs(T),Esd(T)),x0,A,b,Aeq,beq,lb,ub,@(x)mycon2(x,btk(T),stk(T),snum(T),bnum(T),stpv(T),btpv(T),Esd(T)));
    Flag(T,2)=exitflag;
    
    if(btk(T)/price1(1)+stk(T)/price1(2)-snum(T)-bnum(T)-stpv(T)-btpv(T)+Esd(T)>=0)
        pro(T)=pro2;
        price(T,:)=price2;
    elseif(btk(T)/price2(1)+stk(T)/price2(2)-snum(T)-bnum(T)-stpv(T)-btpv(T)+Esd(T)<=0)
        pro(T)=pro1;
        price(T,:)=price1;
    else
        if(pro1<pro2)
            pro(T)=pro1;
            price(T,:)=price1;
        else
            pro(T)=pro2;
            price(T,:)=price2;
        end
    end
    
    for i=1:usernum
        if(role(T,i)==1)
            optload(T,i)=k(T,i)/price(T,1)-1;
            finalprices(T,i)=price(T,1);            
        else
            optload(T,i)=k(T,i)/price(T,2)-1;
            finalprices(T,i)=-price(T,2);            
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
        if(optload(T,i)-pv(T,i)>=0)
            optEsell(T)=optEsell(T)+optload(T,i)-pv(T,i);
        else
            optEbuy(T)=optEbuy(T)-(optload(T,i)-pv(T,i));
        end
    end
end

for T=1:24
    if(-optNLt(T)-Esd(T)>0)
        optMGOpro(T)=price(T,1)*optEsell(T)-price(T,2)*optEbuy(T)+pgb(T)*(-optNLt(T));
    else
        optMGOpro(T)=price(T,1)*optEsell(T)-price(T,2)*optEbuy(T)+pgs(T)*(-optNLt(T));
    end
end
%%

for T=1:24
    for i=1:usernum
        userpro(T,i)=k(T,i)*log(1+l_real(T,i))-pgs(T)*max(l_real(T,i)-pv(T,i),0)-pgb(T)*min(l_real(T,i)-pv(T,i),0)+pvbt*pv(T,i);
        cu(T,i)=k(T,i)*log(1+l_real(T,i));
        tu(T,i)=-pgs(T)*max(l_real(T,i)-pv(T,i),0)-pgb(T)*min(l_real(T,i)-pv(T,i),0)+pvbt*pv(T,i);
        newuserpro(T,i)=k(T,i)*log(1+newload(T,i))-price(T,1)*max(newload(T,i)-pv(T,i),0)-price(T,2)*min(newload(T,i)-pv(T,i),0)+pvbt*pv(T,i);
        newcu(T,i)=k(T,i)*log(1+newload(T,i));
        newtu(T,i)=-price(T,1)*max(newload(T,i)-pv(T,i),0)-price(T,2)*min(newload(T,i)-pv(T,i),0)+pvbt*pv(T,i);
    end
end

for i=1:usernum
    NL(:,i)=l_real(:,i)-pv(:,i);
end

for i=1:usernum
    newNL(:,i)=newload(:,i)-pv(:,i);
end



figure(1)
plot(1:1:24,newuserpro-userpro);
plot(1:1:24,newuserpro-userpro);
pro_agent=sum(optMGOpro);
benefit=sum(newuserpro-userpro,1);%./sum(userpro,1);
benefit_per=100*sum(newuserpro-userpro,1)./sum(userpro,1);