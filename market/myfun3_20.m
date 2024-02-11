function f=myfun3(x,role,stk,btk,kall,bstv,snum,bnum,stpv,btpv,pgs,Esd)

f=stk+btk-sum(x(1,:).*role)-sum(x(1,:).*bstv(1,:))-sum(x(2,:).*(1-role))-sum(x(2,:).*bstv(2,:))+pgs*(stpv+btpv+snum+bnum-sum(kall(2,:)./x(1,:))-sum(kall(3,:)./x(2,:))-Esd);

f=-f;

end