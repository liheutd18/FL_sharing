function f=myfuncvar(x,tl,tpv,wb,ws,rho)

temp=0;
cost=0;

for i=1:125
    temp=temp+max(max(abs(tl(:,i)-tpv(:,i)+x(1:24)))-x(25),0);
    cost=cost+wb*max(tl(:,i)-tpv(:,i)+x(1:24),0)-ws*min(tl(:,i)-tpv(:,i)+x(1:24),0);
end
    f=cost/125+rho*(x(25)+20*temp/125);


% f=wb*max(tl-tpv+x(1:24),0)-ws*min(tl-tpv+x(1:24),0)+rho*(x(25)+20*max(sum(abs(tl-tpv+x(1:24)))-x(25),0));

% f=0;
% for i=1:125
%     f=f+(wb*max(tl(:,i)-tpv(:,i)+x(1:24),0)-ws*min(tl(:,i)-tpv(:,i)+x(1:24),0))/125+rho*20*max(sum(abs(tl(:,i)-tpv(:,i)+x(1:24)))-x(25),0)/125;
% end
% 
% f=f+rho*x(25);
end
