tic;
usernum=3;

pgs=[0.3658 0.3658 0.3658 0.3658 0.3658 0.3658 0.3658 0.8595 0.8595 0.8595 1.3782 1.3782 1.3782 1.3782 1.3782 0.8595 0.8595 0.8595 1.3782 1.3782 1.3782 0.8595 0.8595 0.3658];
pgb=[0.35 0.35 0.35 0.35 0.35 0.35 0.35	0.35 0.35 0.35 0.35	0.35 0.35 0.35 0.35	0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35];
% pgs=1.0*ones(1,24);
% pgb=0.4*ones(1,24);
% perror = -0.05 + 0.1*rand(1,24);
load perror;

a=zeros(24,24);
SOCmin=0.2;
SOCmax=0.9;
SOCini=0.4;
SOCsd=zeros(25,1);
SDcap=5000;
SOCsd(1)=SOCini*SDcap;
c=ones(24,1);

[pv,l,loadmin,loadmax]=ini();

fc=0;
avg=0;
temp=[-1.31472350348688;-0.416411219699434;1.22468782478534;-0.0435842055463326;0.582423277447969;-1.00650007461934;0.0645167423111044;0.600291949185784;-1.36151495486457;0.347592631960065;-0.181843218459334;-0.939534765941492;-0.0375331888191716;-1.89630449362245;-2.12797676818267;-1.17692333071496;-0.990532220484176;-1.17303232726741;-1.72542778952869;0.288228089665011;-1.59418372026681;0.110218849223381;0.787066676357980;-0.00222678631383613];
error=sqrt(fc)*temp+avg;

load pverror;
load lerror;

pvyc=pv+pverror;
lyc=l+lerror;
tl=sum(l,2);
tpv=sum(pv,2);
tlyc=sum(lyc,2);
tpvyc=sum(pvyc,2);
iniNLt=tl-tpv;

for i=1:24
    for j=1:i
        a(i,j)=1;
    end
end
b=-a;
A=[a;b];

x0=zeros(24,1); 
b=[SDcap*(SOCmax-SOCini)*c;SDcap*(SOCini-SOCmin)*c]; Aeq=[]; beq=[];
lb=-2000*c;
ub=2000*c;

wb=ones(1,24);
ws=ones(1,24);
for T=1:24
    wb(T)=pgs(T)/max(pgs);
    ws(T)=min(pgb)/pgb(T);
end

[Esd,cost]=fmincon(@(x)myfuntest(x,tlyc,tpvyc,wb,ws),x0,A,b,Aeq,beq,lb,ub);
% realcost=wb*max(tl-realtpv+Esd,0)-ws*min(tl-realtpv+Esd,0);

toc;
for i=2:25
    SOCsd(i)=SOCsd(i-1)+Esd(i-1);
end

Esd=[Esd;Esd(24)];
SOCsd=SOCsd/SDcap;

figure(1)
plotyy(0:1:24,Esd,0:1:24,SOCsd);