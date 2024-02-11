function f=myfuntest(x,tnl,wb,ws)
eff=1;
%eff=1;
a=length(tnl);
for t=1:a
gamma(t)=1^t;%weak effect
end
f=gamma.*wb*max(tnl+x/eff,0)+gamma.*ws*min(tnl+x*eff,0);
%f=[f;f1];
end