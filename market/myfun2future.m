function [f,r]=myfun2future(newload,l,pv,price,stk,btk,snum,bnum,stpv,btpv,pgb,pgs,Esd)

ab=length(Esd);
dr=l(1,:)-newload;
avedr=dr/ab;
l=l-avedr;
l(1,:)=newload;

ntl=l-pv;
% for T=1:ab
%     ntl(T,:)=ntl(T,:).*((0.1*rand(1,6)-0.05)+1);
% end
x1=price(1)*ones(ab,1);x2=price(2)*ones(ab,1);
x1(2:end)=0.8*x1(2:end);x2(2:end)=1.2*x2(2:end);
% x1=[price(1);0.8*pgs(2:end)'];
% x2=[price(2);1.2*pgb(2:end)'];
% %x1=(pgs)';x2=(pgb)';
% % ntllstm=sum(ntl,2)';
% % ntllstm(2:3,:)=1;
% % load lstmnet0203         
% %         numTimeStepsTest = ab*4;
% %         for i = 2:numTimeStepsTest
% %             [net,ntllstm(:,i)] = predictAndUpdateState(net,ntllstm(:,i-1),'ExecutionEnvironment','cpu');
% %         end
% %         ntllstm = sig.*ntllstm + mu;
% %         ntllstm = ntllstm(:,2:end);
% %DeltaE=stpv+btpv+snum+bnum-(btk./x1)-(stk./x2)-Esd;%%max(DeltaE,0) min(DeltaE,0)
usernum=10;bamount=zeros(ab,1);samount=zeros(ab,1);bnum=zeros(ab,1);snum=zeros(ab,1);
    for T=1:ab
        for i=1:usernum
            if(ntl(T,i)>0)
                role(T,i)=1;
                bnum(T)=bnum(T)+1;
                bamount(T)=bamount(T)+ntl(T,i);
                %bamount(T)=bamount(T)+pv(T,i);
            else
                role(T,i)=0;
                snum(T)=snum(T)+1;
                samount(T)=samount(T)+ntl(T,i);
                %samount(T)=samount(T)+pv(T,i);
            end
        end
        degra(T)=0.95^(T-1);
    end
    
    DeltaE=bamount+samount+Esd;%%max(DeltaE,0) min(DeltaE,0)
%     
%     %f=sum(x1(1)*bamount(1))+sum(x2(1)*samount(1))+(pgs(2:end)*bamount(2:end))+sum(pgb(2:end)*samount(2:end)) -pgs*max(DeltaE,0)-pgb*min(DeltaE,0)-0.15*sum(max(Esd,0));
%     f=x1.*bamount+x2.*samount-pgs'.*max(DeltaE,0)-pgb'.*min(DeltaE,0)-0.15*(max(Esd,0));%-x1.*min(Esd,0)-x2.*max(Esd,0);
     f=x1(1)*bamount(1)+x2(1)*samount(1)-sum(pgs'.*max(DeltaE,0))-sum(pgb'.*min(DeltaE,0))-sum(0.15*(max(Esd,0)));%-x1.*min(Esd,0)-x2.*max(Esd,0);
     
     r=f(1); 
%     f=f(2:end);
     %f=f.*degra';
%     %f=f(1);
     f=sum(f);
%     %f=stk+btk-x1.*(bnum+btpv)-x2.*(snum+stpv)+pgb*max(DeltaE,0)+pgs*min(DeltaE,0);
%     quantity=[bamount,samount,Esd];
    %nowf=price(1)*quantity(1,1)+price(2)*quantity(1,2)-pgs(T)*max(quantity(1,3),0)-pgb(T)*min(quantity(1,3),0);
    %f=f;%+nowf;
%f=sum(-f);

end