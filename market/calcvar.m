usernum=3;

pgs=[0.3658 0.3658 0.3658 0.3658 0.3658 0.3658 0.3658 0.8595 0.8595 0.8595 1.3782 1.3782 1.3782 1.3782 1.3782 0.8595 0.8595 0.8595 1.3782 1.3782 1.3782 0.8595 0.8595 0.3658];
pgb=[0.35 0.35 0.35 0.35 0.35 0.35 0.35	0.35 0.35 0.35 0.35	0.35 0.35 0.35 0.35	0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35];
% pgs=1.0*ones(1,24);
% pgb=0.4*ones(1,24);

a=zeros(24,24);
SOCmin=0.2;
SOCmax=0.9;
SOCini=0.4;
SOCsd=zeros(25,1);
SDcap=5000;
SOCsd(1)=SOCini*SDcap;
c=ones(24,1);

[pv,l,lmin,lmax]=ini();

% fc=0;
% avg=0;
% temp=[-1.31472350348688;-0.416411219699434;1.22468782478534;-0.0435842055463326;0.582423277447969;-1.00650007461934;0.0645167423111044;0.600291949185784;-1.36151495486457;0.347592631960065;-0.181843218459334;-0.939534765941492;-0.0375331888191716;-1.89630449362245;-2.12797676818267;-1.17692333071496;-0.990532220484176;-1.17303232726741;-1.72542778952869;0.288228089665011;-1.59418372026681;0.110218849223381;0.787066676357980;-0.00222678631383613];
% error=sqrt(fc)*temp+avg;

tl=sum(l,2);
tpv=sum(pv,2);
% realtpv=tpv+error;

% r1 = -1000 + 2000*rand(24,125);
% r2 = -1000 + 2000*rand(24,125);
load r1;
load r2;
for i=1:125
    rl(:,i)=r1(:,i)+tl;
    rpv(:,i)=r2(:,i)+tpv;
end


for i=1:24
    for j=1:i
        a(i,j)=1;
    end
end
b=-a;
A=[a;b];
yzero=zeros(48,1);
A=[A,yzero];
xzero=zeros(1,25);
A=[A;xzero];

x0=zeros(25,1); 
b=[SDcap*(SOCmax-SOCini)*c;SDcap*(SOCini-SOCmin)*c];
b=[b;0];
Aeq=[]; beq=[];
lb=-2000*c;
ub=2000*c;
lb=[lb;-inf];
ub=[ub;inf];

wb=ones(1,24);
ws=ones(1,24);
for T=1:24
    wb(T)=pgs(T)/max(pgs);
    ws(T)=min(pgb)/pgb(T);
end

rho=1;

[Esd,fvalue]=fmincon(@(x)myfuncvar(x,rl,rpv,wb,ws,rho),x0,A,b,Aeq,beq,lb,ub);
temp=0;
cost=0;
for i=1:125
    temp=temp+max(max(abs(rl(:,i)-rpv(:,i)+Esd(1:24)))-Esd(25),0);
    cost=cost+wb*max(rl(:,i)-rpv(:,i)+Esd(1:24),0)-ws*min(rl(:,i)-rpv(:,i)+Esd(1:24),0);
end
cost=cost/125;
cvar=Esd(25)+20*temp/125;
%realcost=wb*max(tl-realtpv+Esd,0)-ws*min(tl-realtpv+Esd,0);

% for i=2:25
%     SOCsd(i)=SOCsd(i-1)+Esd(i-1);
% end
% 
% Esd=[Esd;Esd(24)];
% SOCsd=SOCsd/SDcap;
% 
% figure(1)
% plotyy(0:1:24,Esd,0:1:24,SOCsd);
nl=zeros(24,125);
for i=1:125
    nl(:,i)=rl(:,i)-rpv(:,i)+Esd(1:24);
end
maxnl=max(nl);
theone=max(maxnl);