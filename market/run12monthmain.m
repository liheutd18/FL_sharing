for capacitym=1:50
%for startday=1:365

% clc 
% clear
% close all
usernum=10;

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
%capacitym=40;
SDcap=capacitym*1000;
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

%[pv,l,loadmin,loadmax]=ini();
[pv_real,l_real,loadmin,loadmax,estimatepv]=Pecan();
startday=313; %31 310  151
pv=pv_real((startday-1)*24+1:24*startday,:);l=l_real((startday-1)*24+1:24*startday,:);
%      Data=readtable('Forecasting results.xlsx','Sheet','sheet3');
%      for j=1:10
%          NetLoad_real(j,:)=table2array(Data(2+(j-1)*4,50:73)); 
%          NetLoad_pred(j,:)=table2array(Data(1+(j-1)*4,50:73)); 
%      end
%          AGGNetLoad_predDA=table2array(Data(42,50:73)); 
%          AGGNetLoad_real=table2array(Data(43,50:73));      
%          AGGNetLoad_predHA=table2array(Data(46,50:73)); 
%      pv=estimatepv((30-1)*24+1:30*24,:);
%      pv_real=pv_real((30-1)*24+1:30*24,:);% day =30
pv(pv<0)=0;pv=pv*1000;
%pv_real(pv_real<0)=0;pv_real=pv_real*1000;

 l=l;l=l*1000;
% loadmin=loadmin*1000;loadmax=loadmax*1000;
% ntl=NetLoad_real';ntl=ntl*1000;
% ntl_pred=NetLoad_pred';ntl_pred=ntl_pred*1000;
% 
% l=ntl+pv;

% l=ntl_pred+pv;
% l(l<=200)=200;
% l_real=ntl+pv_real;
loadmin=0.8*l;
loadmax=1.2*l;
% tl=sum(l,2);
% tpv=sum(pv,2);
%iniNLt=tl-tpv;

%tnl=AGGNetLoad_predDA';tnl=tnl*1000;
tl=sum(l,2);
tpv=sum(pv,2);
tnl=sum(l-pv,2);
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
%tic;
%     wb=[0.2*ones(1,7) 0.55*ones(1,16)  0.2*ones(1,1) ];
%     ws=0.1*ones(1,24);
    c_rate=2;%1%2
%    for i=1:500
%        SDcap=40*1000;
        b=[SDcap*(SOCmax-SOCini)*c;SDcap*(SOCini-SOCmin)*c];
        lb=-(SDcap/c_rate)*c;
        ub=(SDcap/c_rate)*c;
        [Esd,cost,cost1]=fmincon(@(x)myfuntest(x,tnl,wb,ws),x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
%         cost=cost1+0.15*i;
%    end
%toc;
New_NL=max(tnl+Esd/eff,0)+min(tnl+Esd*eff,0);

for i=2:25
    SOCsd(i)=SOCsd(i-1)+Esd(i-1);
end
SOCsdall(:,startday)=SOCsd;
%end
%      K_means= kmeans(SOCsdall',4);[row, col] = find(isnan(K_means));
%         l=1;m=1;n=1;o=1;
%         for i=1:365
%             if K_means(i)==1
%                K_1(l,:) =SOCsdall(:,i);
%                l=l+1;
%             elseif K_means(i)==2
%                K_2(m,:) =SOCsdall(:,i);
%                m=m+1;
%             elseif K_means(i)==3
%                K_3(n,:) =SOCsdall(:,i);
%                n=n+1;
%             elseif K_means(i)==4
%                K_4(o,:) =SOCsdall(:,i);
%                o=o+1;
%             end
%         end
%         figure
% %         plot(K_1')
% %         hold on
%         plot(mean(K_1,1),'b--o')
%         hold on
% %         plot(K_2')
% %         hold on
%         plot(mean(K_2,1),'r--o')
%         hold on
% %         plot(K_3')
% %         hold on
%         plot(mean(K_3,1),'g--o');        
%         hold on
% %         plot(K_4')
% %         hold on
%         plot(mean(K_4,1),'y--o');

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
    %%k
%     for i=1:usernum
%         if(l_real(T,i)>=pv(T,i))
%             k(T,i)=(l_real(T,i)+1)*pgs(T);
%              role(T,i)=1;
% %             bnum(T)=bnum(T)+1;
% %             btpv(T)=btpv(T)+pv(T,i);
% %             btk(T)=btk(T)+k(T,i);
% %             bk(T,i)=k(T,i);
% %             btv(T,i)=pv(T,i);
%         else
%             k(T,i)=(l_real(T,i)+1)*pgb(T);
%              role(T,i)=0;
% %             snum(T)=snum(T)+1;
% %             stpv(T)=stpv(T)+pv(T,i);
% %             stk(T)=stk(T)+k(T,i);
% %             sk(T,i)=k(T,i);
% %             stv(T,i)=pv(T,i);
%         end
%     end
    
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
        if(role(T,i)==1)
            optload(T,i)=k(T,i)/priceb(T,i)-1;
            finalprices(T,i)=priceb(T,i);
        else
            optload(T,i)=k(T,i)/prices(T,i)-1;
            finalprices(T,i)=-prices(T,i);
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

for T=1:24
    if(+newNLt(T)+Esd(T)<0)
        %optMGOpro(T)=price(T,1)*optEsell(T)-price(T,2)*optEbuy(T)+pgb(T)*(-optNLt(T));
        optMGOpro(T)=sum((optload(T,:)-pv(T,:)).*finalprices(T,:))+pgb(T)*(-optNLt(T)-Esd(T));
    else
        %optMGOpro(T)=price(T,1)*optEsell(T)-price(T,2)*optEbuy(T)+pgs(T)*(-optNLt(T));
        optMGOpro(T)=sum((optload(T,:)-pv(T,:)).*finalprices(T,:))+pgs(T)*(-optNLt(T)-Esd(T));
    end
end


for T=1:24
    for i=1:usernum
        userpro(T,i)=k(T,i)*log(1+l(T,i))-pgs(T)*max(l(T,i)-pv(T,i),0)-pgb(T)*min(l(T,i)-pv(T,i),0)+pvbt*pv(T,i);
        cu(T,i)=k(T,i)*log(1+l(T,i));
        tu(T,i)=-pgs(T)*max(l(T,i)-pv(T,i),0)-pgb(T)*min(l(T,i)-pv(T,i),0)+pvbt*pv(T,i);
        
        newuserpro(T,i)=k(T,i)*log(1+newload(T,i))-finalprices(T,i)*max(newload(T,i)-pv(T,i),0)+finalprices(T,i)*min(newload(T,i)-pv(T,i),0)+pvbt*pv(T,i);
%         newcu(T,i)=k(T,i)*log(1+newload(T,i));
%         newtu(T,i)=-finalprices(T,i)*max(newload(T,i)-pv(T,i),0)-finalprices(T,2)*min(newload(T,i)-pv(T,i),0)+pvbt*pv(T,i);
    end
end

% for i=1:usernum
%     NL(:,i)=l(:,i)-pv(:,i);
% end
% 
% for i=1:usernum
%     newNL(:,i)=newload(:,i)-pv(:,i);
% end


plot(priceb)
hold on
plot(prices)
plot(finalprices)
plot(1:1:24,newuserpro-userpro);
pro_agent=sum(optMGOpro)/10000;

benefit=sum(newuserpro-userpro,1);%./sum(userpro,1);
benefit_per=100*sum(newuserpro-userpro,1)./sum(userpro,1);
benefit_per1(capacitym,:)=100*sum(newuserpro-userpro,1)./sum(userpro,1);
uti1=sum(newuserpro,1)/10000;
uti2=sum(cu+tu,1)/10000;
display(pro_agent)

end